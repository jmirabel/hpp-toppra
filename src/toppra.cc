// Copyright (c) 2020, Joseph Mirabel
// Authors: Joseph Mirabel (joseph.mirabel@laas.fr)
//          Olivier Roussel (olivier.roussel@laas.fr)
//

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.

#include "toppra.hh"

#include <sstream>

#include <hpp/core/path-optimizer.hh>
#include <hpp/core/path-vector.hh>
#include <hpp/core/path.hh>
#include <hpp/core/plugin.hh>
#include <hpp/core/problem-solver.hh>
#include <hpp/core/problem.hh>
#include <hpp/core/time-parameterization/piecewise-polynomial.hh>
#include <hpp/pinocchio/device.hh>
#include <hpp/pinocchio/joint-collection.hh>
#include <pinocchio/multibody/model.hpp>
#include <toppra/algorithm/toppra.hpp>
#include <toppra/constraint/joint_torque/pinocchio.hpp>
#include <toppra/constraint/linear_joint_acceleration.hpp>
#include <toppra/constraint/linear_joint_velocity.hpp>
#include <toppra/geometric_path.hpp>
#include <toppra/solver/glpk-wrapper.hpp>
#include <toppra/solver/qpOASES-wrapper.hpp>
#include <toppra/solver/seidel.hpp>
#include <toppra/toppra.hpp>

#include "piecewise-polynomial.hh"
#include "serialization.hh"

namespace hpp {
namespace core {

class PathWrapper : public toppra::GeometricPath {
 public:
  PathWrapper(PathPtr_t path)
      : toppra::GeometricPath((int)path->outputSize(),
                              (int)path->outputDerivativeSize()),
        path_(path) {}

  toppra::Vector eval_single(toppra::value_type time, int order) const {
    bool success;
    toppra::Vector res;
    if (order == 0) {
      res = path_->eval(time, success);
      assert(success);
    } else {
      res.resize(dof());
      path_->derivative(res, time, order);
    }
    return res;
  }

  toppra::Bound pathInterval() const {
    const interval_t& tr = path_->timeRange();
    return (toppra::Bound() << tr.first, tr.second).finished();
  }

 private:
  PathPtr_t path_;
};

namespace timeParameterization {
class Extract : public TimeParameterization {
 public:
  Extract(TimeParameterizationPtr_t inner, value_type dt, value_type ds)
      : inner_(inner), dt_(dt), ds_(ds) {}

  value_type value(const value_type& t) const {
    return inner_->value(t + dt_) + ds_;
  }

  value_type derivative(const value_type& t, const size_type& order) const {
    return inner_->derivative(t + dt_, order);
  }

  value_type derivativeBound(const value_type& low,
                             const value_type& up) const {
    return derivativeBound(low + dt_, up + dt_);
  }

  TimeParameterizationPtr_t copy() const {
    return TimeParameterizationPtr_t(new Extract(inner_, dt_, ds_));
  }

 private:
  TimeParameterizationPtr_t inner_;
  value_type dt_, ds_;
};
}  // namespace timeParameterization

#define PARAM_HEAD "PathOptimization/TOPPRA/"

namespace pathOptimization {

TimeParameterizationPtr_t constantAccelerationParametrization(
    toppra::Vector const& t,
    toppra::Vector const& s,
    toppra::Vector const& sd)
{
  // Inputs
  // int N
  // s, sd, t
  const size_type N = t.size() - 1;

  // P_i(t) = c_0 + c_1 t + c_2 t**2
  // P_i(t[i]) = s[i]
  // P_i(t[i+1]) = s[i+1]
  // P_i'(t[i]) = sd[i]
  // P_i'(t[i+1]) = sd[i+1]
  // P_i''(.) = u[i] = (sd[i+1]**2 - sd[i]**2) / (2 * (s[i+1]-s[i]))
  //
  // c_2 = 0.5 u[i]
  // c_1 + 2 c_2 t[i] = sd[i]
  // c_0 + c_1 t[i] + c_2 t[i]**2 = s[i]
  constexpr int order = 2;
  typedef timeParameterization::ShiftedPiecewisePolynomial<order> TimeParam_t;
  auto c = TimeParam_t::ParameterMatrix_t(order + 1, N);
  for (size_type i = 0; i < N; ++i) {
    c(2, i) = (sd[i+1] - sd[i]) * (sd[i+1] + sd[i]) / (s[i+1] - s[i]) / 4;
    c(1, i) = sd[i];
    c(0, i) = s[i];
  }

  return std::make_shared<TimeParam_t>(c, t);
}

TimeParameterizationPtr_t hermiteCubicSplineParametrization(
    toppra::Vector const& t,
    toppra::Vector const& s,
    toppra::Vector const& sd)
{
  // Inputs
  // int N
  // s, sd, t

  const size_type N = t.size() - 1;

  constexpr int order = 3;
  typedef timeParameterization::PiecewisePolynomial<order> TimeParam_t;
  auto spline_coeffs = TimeParam_t::ParameterMatrix_t(order + 1, N);
  for (size_type i = 1; i <= N; ++i) {
    const auto inv_dt = 1. / (t[i] - t[i - 1]);
    const auto inv_dt2 = inv_dt * inv_dt;
    const auto inv_dt3 = inv_dt2 * inv_dt;
    const auto t_p = t[i - 1];
    const auto t_p2 = t_p * t_p;
    const auto t_p3 = t_p2 * t_p;
    const auto ds = (s[i] - s[i - 1]);
    const auto b = (2 * sd[i - 1] + sd[i]) * inv_dt;
    const auto c = (sd[i - 1] + sd[i]) * inv_dt2;
    spline_coeffs(0, i - 1) = 2 * t_p3 * ds * inv_dt3 +
                              3 * t_p2 * ds * inv_dt2 - t_p3 * c - t_p2 * b -
                              t_p * sd[i - 1] + s[i - 1];
    spline_coeffs(1, i - 1) = -6 * t_p2 * ds * inv_dt3 -
                              6 * t_p * ds * inv_dt2 + 3 * t_p2 * c +
                              2 * t_p * b + sd[i - 1];
    spline_coeffs(2, i - 1) =
        6 * t_p * ds * inv_dt3 + 3 * ds * inv_dt2 - 3 * t_p * c - b;
    spline_coeffs(3, i - 1) = -2 * ds * inv_dt3 + c;
  }

  return std::make_shared<TimeParam_t>(spline_coeffs, t);
}

/// \param[out] id_subpaths \c id_subpaths[i] the index in gridpoints that
///             corresponds to the start of \c paths[i]
toppra::Vector evenlySpacedGridpoints(
    std::vector<PathPtr_t> const& paths,
    size_type const N,
    value_type const maxSegmentLength,
    std::vector<size_type>& id_subpaths)
{
  std::vector<value_type> S;
  S.reserve(N + paths.size() - 1);

  id_subpaths.reserve(paths.size() + 1);

  S.push_back(0.);
  id_subpaths.push_back(0);
  for (PathPtr_t const& subpath : paths) {
    auto I = subpath->paramRange();
    value_type pathS = I.second - I.first;
    size_type n = size_type(std::ceil(pathS / maxSegmentLength));
    value_type p0 = S.back();
    for (size_type k = 1; k <= n; ++k)
      S.push_back(p0 + (pathS * (value_type)k) / (value_type)n);
    id_subpaths.push_back(S.size() - 1);
  }
  return toppra::Vector(Eigen::Map<toppra::Vector>(S.data(), S.size()));
}

/// \param[out] id_subpaths \c id_subpaths[i] the index in gridpoints that
///             corresponds to the start of \c paths[i]
toppra::Vector evenlyTimeSpacedGridpoints(
    std::shared_ptr<PathWrapper> pathWrapper,
    std::vector<PathPtr_t> const& paths,
    size_type const N,
    value_type const maxSegmentLength,
    std::vector<size_type>& id_subpaths)
{
  toppra::Vector initialS(paths.size() + 1);

  //id_subpaths.reserve(paths.size() + 1);

  initialS[0] = 0.0;
  //id_subpaths.push_back(0);
  for (auto i = 0ul; i < paths.size(); ++i)
    initialS[i+1] = initialS[i] + paths[i]->length();

  const double maxErrorThreshold = 1e-4;
  const int maxIterations = 100;
  toppra::Vector gridpoints (pathWrapper->proposeGridpoints(
      maxErrorThreshold,
      maxIterations,
      maxSegmentLength,
      static_cast<int>(N),
      initialS));
  // id_subpaths[i] is the index in gridpoints that corresponds to the start of paths[i], aka initialS[i]
  int k = 0;
  for (int i = 0; i < initialS.size(); ++i) {
    bool found = false;
    for (; k < gridpoints.size(); ++k) {
      if (gridpoints[k] == initialS[i]) {
        found = true;
        break;
      }
    }
    if (!found)
      throw std::logic_error("Initial gridpoints not found in proposed gridpoints.");
    id_subpaths[i] = k;
  }
  return gridpoints;
}

void TOPPRA::inputSerialization(PathPtr_t path) const {
  const std::string filename =
      problem()->getParameter(PARAM_HEAD "inputSerialization").stringValue();
  if (filename.size() == 0)
    return;

  DevicePtr_t device(problem()->robot());
  if (filename.substr(filename.size() - 4) == ".txt")
    parser::serializePath<serialization::text_oarchive>(device, path, filename);
  else
    parser::serializePath<serialization::binary_oarchive>(device, path, filename);
}

toppra::LinearConstraintPtrs TOPPRA::constraints()
{
  const value_type effortScale =
      problem()->getParameter(PARAM_HEAD "effortScale").floatValue();
  const value_type velScale =
      problem()->getParameter(PARAM_HEAD "velocityScale").floatValue();
  const vector_t accLimits =
      problem()->getParameter(PARAM_HEAD "accelerationLimits").vectorValue();

  const pinocchio::Model& model = problem()->robot()->model();

  using namespace toppra::constraint;

  // Create the TOPPRA constraints
  toppra::LinearConstraintPtrs v;

  // Joint velocity limits
  v.push_back(
      std::make_shared<LinearJointVelocity>(-velScale * model.velocityLimit,
                                            velScale * model.velocityLimit)
  );
  // Joint acceleration limits
  if (accLimits.size() > 0) {
    if (accLimits.size() != model.nv) {
      std::ostringstream oss;
      oss << "Acceleration limits should be of size " << model.nv << " and a "
        "vector of size " << accLimits.size() << " is provided.";
      throw std::invalid_argument(oss.str());
    }
    v.push_back(std::make_shared<LinearJointAcceleration>(-accLimits, accLimits));
  }
  // Joint torque limits
  if (effortScale >= 0) {
    auto torqueConstraint =
        std::make_shared<jointTorque::Pinocchio<pinocchio::Model> >(model);  // No friction
    torqueConstraint->lowerBounds(effortScale * torqueConstraint->lowerBounds());
    torqueConstraint->upperBounds(effortScale * torqueConstraint->upperBounds());
    v.push_back(torqueConstraint);
  }
  for (auto& c : v)
    c->discretizationType(toppra::Interpolation);

  return v;
}

TOPPRA::InterpolationMethod TOPPRA::interpolationMethod() const
{
  const std::string interpolationMethod =
      problem()->getParameter(PARAM_HEAD "interpolationMethod").stringValue();
  if (interpolationMethod == "hermite") {
    return Hermite;
  } else if (interpolationMethod == "constant_acceleration") {
    return ConstantAcceleration;
  } else {
    std::ostringstream oss;
    oss << "Invalid interpolationMethod. Allowed values are 'hermite' and "
      "'constant_acceleration'. Provided value: " << interpolationMethod;
    throw std::invalid_argument(oss.str());
  }
}

TOPPRA::GridpointMethod TOPPRA::gridpointMethod() const
{
  const std::string gridpointMethod =
      problem()->getParameter(PARAM_HEAD "gridpointMethod").stringValue();
  if (gridpointMethod == "param_space") {
    return EvenlyParamSpaced;
  } else if (gridpointMethod == "time_space") {
    return EvenlyTimeSpaced;
  } else {
    std::ostringstream oss;
    oss << "Invalid gridpointMethod (method to generate gridpoints). Allowed values are 'param_space' and "
      "'time_space'. Provided value: " << gridpointMethod;
    throw std::invalid_argument(oss.str());
  }
}

PathVectorPtr_t TOPPRA::optimize(const PathVectorPtr_t& path)
{
  inputSerialization(path);

  const size_type solver =
      problem()->getParameter(PARAM_HEAD "solver").intValue();

  toppra::LinearConstraintPtrs v = std::move(constraints());

  PathVectorPtr_t flatten_path =
      PathVector::create(path->outputSize(), path->outputDerivativeSize());
  path->flatten(flatten_path);

  const value_type min_path_length = 1e-6;
  if (path->length() < min_path_length) {
    return flatten_path;
  }

  size_type N = problem()->getParameter(PARAM_HEAD "N").intValue();

  std::vector<PathPtr_t> paths(flatten_path->numberPaths());
  for (auto i = 0ul; i < flatten_path->numberPaths(); ++i)
    paths[i] = flatten_path->pathAtRank(i);
  value_type maxSegmentLength = flatten_path->length() / (value_type)N;

  std::shared_ptr<PathWrapper> pathWrapper(new PathWrapper(flatten_path));

  // 1. Compute TOPPRA grid points (in the parameter space).
  std::vector<size_type> id_subpaths;
  id_subpaths.reserve(flatten_path->numberPaths() + 1);

  toppra::Vector gridpoints;
  switch(gridpointMethod()) {
    case EvenlyTimeSpaced:
      gridpoints = std::move(evenlyTimeSpacedGridpoints(
            pathWrapper, paths, N, maxSegmentLength, id_subpaths));
      break;
    case EvenlyParamSpaced:
      gridpoints = std::move(evenlySpacedGridpoints(paths, N, maxSegmentLength, id_subpaths));
      break;
  }
  N = gridpoints.size() - 1;

  // 2. Apply TOPPRA on the full path
  toppra::algorithm::TOPPRA algo(v, pathWrapper);
  algo.setN((int)N);
  switch (solver) {
    default:
      hppDout(error, "Solver " << solver << " does not exists. Using Seidel");
    case 0:
      algo.solver(std::make_shared<toppra::solver::Seidel>());
      break;
    case 1:
      algo.solver(std::make_shared<toppra::solver::GLPKWrapper>());
      break;
    case 2:
      algo.solver(std::make_shared<toppra::solver::qpOASESWrapper>());
      break;
  }
  algo.setGridpoints(gridpoints);
  auto ret_code = algo.computePathParametrization();
  if (ret_code != toppra::ReturnCode::OK) {
    std::stringstream ss;
    ss << "TOPPRA failed, returned code: " << static_cast<int>(ret_code)
       << std::endl;
    throw std::runtime_error(ss.str());
  }
  const auto out_data = algo.getParameterizationData();
  assert(out_data.gridpoints.size() == out_data.parametrization.size());
  toppra::Vector sd(out_data.parametrization.cwiseMax(0.).cwiseSqrt());

  // 3. Build the parameterization function and apply it on each subpaths.
  PathVectorPtr_t res =
      PathVector::create(path->outputSize(), path->outputDerivativeSize());

  // 3.1 forward integration of time parameterization (trapezoidal integration)
  vector_t t = vector_t(N + 1);
  t[0] = 0.;  // start time is 0
  const auto& s = out_data.gridpoints;
  for (size_type i = 1; i <= N; ++i) {
    const auto sd_avg = (sd[i - 1] + sd[i]) * 0.5;
    const auto ds = s[i] - s[i - 1];
    assert(sd_avg > 0.);
    const auto dt = ds / sd_avg;
    t[i] = t[i - 1] + dt;
  }

  // 3.2 time parameterization based on
  // - piecewise constant acceleration or
  // - hermite cubic spline interpolation
  TimeParameterizationPtr_t global;
  switch(interpolationMethod()) {
    case ConstantAcceleration:
      global = constantAccelerationParametrization(t, s, sd);
      break;
    case Hermite:
      global = hermiteCubicSplineParametrization(t, s, sd);
      break;
  }

  for (auto i = 0ul; i < paths.size(); ++i) {
    value_type t0 = t[id_subpaths[i]];
    value_type p0 = -gridpoints[id_subpaths[i]];
    paths[i]->timeParameterization(
        TimeParameterizationPtr_t(
            new timeParameterization::Extract(global, t0, p0)),
        interval_t(0, t[id_subpaths[i + 1]] - t0));
    res->appendPath(paths[i]);
  }

  lastTimeParameterization_ = global;
  return res;
}

HPP_START_PARAMETER_DECLARATION(TOPPRA)
Problem::declareParameter(ParameterDescription(Parameter::VECTOR,
                                               PARAM_HEAD "accelerationLimits",
                                               "Define the acceleration limits.",
                                               Parameter(vector_t())));
Problem::declareParameter(ParameterDescription(Parameter::STRING,
                                               PARAM_HEAD "gridpointMethod",
                                               "Define the method for generating gridpoints.\n"
                                               "Accepted values are:\n"
                                               "  \"param_space\": Evenly spaced in parameter,\n"
                                               "  \"time_space\": Evenly spaced in time",
                                               Parameter(std::string("param_space"))));
Problem::declareParameter(ParameterDescription(Parameter::STRING,
                                               PARAM_HEAD "interpolationMethod",
                                               "Define the interpolation method for the output of TOPPRA.\n"
                                               "Accepted values are: \"hermite\", \"constant_acceleration\"",
                                               Parameter(std::string("constant_acceleration"))));
Problem::declareParameter(ParameterDescription(Parameter::FLOAT,
                                               PARAM_HEAD "effortScale",
                                               "Effort rescaling value.",
                                               Parameter((value_type)-1)));
Problem::declareParameter(ParameterDescription(Parameter::FLOAT,
                                               PARAM_HEAD "velocityScale",
                                               "Velocity rescaling value.",
                                               Parameter((value_type)1)));
Problem::declareParameter(ParameterDescription(Parameter::INT,
                                               PARAM_HEAD "solver",
                                               "0: Seidel\n"
                                               "1: GLPK\n"
                                               "2: qpOASES",
                                               Parameter((size_type)0)));
Problem::declareParameter(ParameterDescription(Parameter::INT, PARAM_HEAD "N",
                                               "Minimal number of sampling point.",
                                               Parameter((size_type)50)));
Problem::declareParameter(ParameterDescription(Parameter::STRING,
                                               PARAM_HEAD "inputSerialization",
                                               "Filename where to serialize the input path. Leave empty to skip.\n"
                                               "Text serialization is used if filename ends with '.txt'.",
                                               Parameter(std::string(""))));
HPP_END_PARAMETER_DECLARATION(TOPPRA)
}  // namespace pathOptimization

class TOPPRAPlugin : public ProblemSolverPlugin {
 public:
  TOPPRAPlugin() : ProblemSolverPlugin("TOPPRAPlugin", "0.0") {}

 protected:
  virtual bool impl_initialize(ProblemSolverPtr_t ps) {
    ps->pathOptimizers.add("TOPPRA", pathOptimization::TOPPRA::create);
    return true;
  }
};
}  // namespace core
}  // namespace hpp

HPP_CORE_DEFINE_PLUGIN(hpp::core::TOPPRAPlugin)
