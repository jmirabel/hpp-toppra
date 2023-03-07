#define CATCH_CONFIG_MAIN
#include"tests.hh"

#include<iostream>
#include<pinocchio/multibody/model.hpp>
#include<hpp/pinocchio/urdf/util.hh>

void print_path_evaluation(hc::PathPtr_t p, int N) {
  hp::vector_t q (p->outputSize()), v(p->outputDerivativeSize()), a(p->outputDerivativeSize());
  hp::value_type
    t = 0.0,
    t0 = p->timeRange().first,
    t1 = p->timeRange().second;
  std::cout
    << "Time interval: " << t0 << ", " << t1
    << "\nt, q, v, a\n";
  for (int i = 0; i <= N; ++i) {
    t = t0 + hp::value_type(i) / N * (t1 - t0);
    p->eval(q, t);
    p->derivative(v, t, 1);
    p->derivative(a, t, 2);
    std::cout << t << ", " << q << ", " << v << ", " << a << '\n';
  }
}

void check_path_velocity(hc::PathPtr_t p) {
  hp::vector_t v(p->outputDerivativeSize());
  hp::value_type
    t = 0.0,
    t0 = p->timeRange().first,
    t1 = p->timeRange().second;
  int N = 1000;
  for (int i = 0; i <= N; ++i) {
    t = t0 + hp::value_type(i) / N * (t1 - t0);
    p->derivative(v, t, 1);
    if (std::abs(v[0]) < 0.01) {
      WARN("small velocity at " << t << ", " << v);
    }
  }
}

static const hc::value_type velLimit = 1.;
static const hc::value_type accLimit = 2.;
static const hc::value_type constraintRelativeTol = 1.02;

hp::DevicePtr_t makeDevice() {
  hp::DevicePtr_t device = hp::Device::create("test");
  hp::urdf::loadModelFromString(device, 0, "", "prismatic_x",
      "<robot name='n'><link name='base_link'/></robot>", "");

  device->model().velocityLimit << velLimit;

  return device;
}

hc::PathVectorPtr_t makeCubicSpline(hc::ProblemPtr_t p) {
  auto sm = SteerSplineOrder3::create(p);
  hp::vector_t q1(p->robot()->neutralConfiguration()), q2(q1);
  hp::matrix_t d1(p->robot()->numberDof(), 1),
               d2(p->robot()->numberDof(), 1);
  std::vector<int> orders {1};

  q1 << 0.0;
  q2 << 1.0;
  d1 << 0.1;
  d2 << 0.1;

  auto spline = sm->steer(q1, orders, d1, q2, orders, d2);
  auto inputPath = hc::PathVector::create(p->robot()->configSize(), p->robot()->numberDof());
  inputPath->appendPath(spline);
  return inputPath;
}

hc::PathVectorPtr_t makeQuinticSpline(hc::ProblemPtr_t p, hc::value_type v0, hc::value_type v1) {
  auto sm = SteerSplineOrder5::create(p);
  hp::vector_t q1(p->robot()->neutralConfiguration()), q2(q1);
  hp::matrix_t d1(p->robot()->numberDof(), 2),
               d2(p->robot()->numberDof(), 2);
  std::vector<int> orders {1, 2};

  q1 << 0.0;
  q2 << 1.0;
  d1 << v0, 0.0;
  d2 << v1, 0.0;

  auto spline = sm->steer(q1, orders, d1, q2, orders, d2);
  auto inputPath = hc::PathVector::create(p->robot()->configSize(), p->robot()->numberDof());
  inputPath->appendPath(spline);
  return inputPath;
}

TEST_CASE("Test throw on invalid parameters") {
  hp::DevicePtr_t device = makeDevice();

  std::string gridpointMethod = GENERATE("time_space", "param_space", "throw: invalid");
  std::string interpolationMethod = GENERATE("hermite", "constant_acceleration", "throw: invalid");
  bool expect_throw = (
      interpolationMethod.compare(0, 5, "throw") == 0
      ||
      gridpointMethod.compare(0, 5, "throw") == 0
  );

  CAPTURE(gridpointMethod, expect_throw);

  auto problem = hc::Problem::create(device);
  problem->setParameter("PathOptimization/TOPPRA/gridpointMethod", hc::Parameter(gridpointMethod));
  problem->setParameter("PathOptimization/TOPPRA/interpolationMethod", hc::Parameter(interpolationMethod));
  hc::PathVectorPtr_t spline = makeCubicSpline(problem);
  auto opt = hc::pathOptimization::TOPPRA::create(problem);
  if (expect_throw) {
    CHECK_THROWS(opt->optimize(spline));
  } else {
    CHECK_NOTHROW(opt->optimize(spline));
  }
}

TEST_CASE("Test parameters") {
  hp::DevicePtr_t device = makeDevice();

  auto problem = hc::Problem::create(device);
  hc::PathVectorPtr_t spline = makeCubicSpline(problem);
  auto opt = hc::pathOptimization::TOPPRA::create(problem);
  hp::vector_t accLimits;

  problem->setParameter("PathOptimization/TOPPRA/interpolationMethod", hc::Parameter(std::string("hermite")));
  CHECK_NOTHROW(opt->optimize(spline));

  problem->setParameter("PathOptimization/TOPPRA/N", hc::Parameter(hp::size_type(200)));
  problem->setParameter("PathOptimization/TOPPRA/interpolationMethod", hc::Parameter(std::string("constant_acceleration")));
  CHECK_NOTHROW(opt->optimize(spline));

  accLimits.resize(2);
  problem->setParameter("PathOptimization/TOPPRA/accelerationLimits", hc::Parameter(accLimits));
  CHECK_THROWS(opt->optimize(spline));

  accLimits.resize(1);
  accLimits << accLimit;
  problem->setParameter("PathOptimization/TOPPRA/accelerationLimits", hc::Parameter(accLimits));
  CHECK_NOTHROW(opt->optimize(spline));
}

void problemParameterToCatch2Info(hc::ProblemPtr_t problem) {
#define PARAM(name, type) \
  INFO(#name ": " << problem->getParameter("PathOptimization/TOPPRA/" #name).type ## Value());

  PARAM(gridpointMethod, string);
  PARAM(interpolationMethod, string);
  PARAM(accelerationLimits, vector);
  PARAM(N, int);

#undef PARAM
}

TEST_CASE("Run TOPPRA optimizer: all cases") {
  hp::DevicePtr_t device = makeDevice();
  auto problem = hc::Problem::create(device);

  std::vector<hc::PathVectorPtr_t> inputPaths = {
      makeCubicSpline(problem),
      makeQuinticSpline(problem, 0.2, 0.2),
      makeQuinticSpline(problem, 0.0, 0.1),
  };

  int iInputPath = GENERATE(0, 1, 2);
  hc::PathVectorPtr_t inputPath = inputPaths[iInputPath];

  std::vector<hc::vector_t> accLimitss = {
    hc::vector_t(),
    (hc::vector_t(1) << accLimit).finished(),
  };
  int iAccLimits = GENERATE(0, 1);
  hc::vector_t accLimits = accLimitss[iAccLimits];

  std::string interpolationMethod = GENERATE("hermite", "constant_acceleration");
  std::string gridpointMethod = GENERATE("time_space", "param_space");

  problem->setParameter("PathOptimization/TOPPRA/accelerationLimits", hc::Parameter(accLimits));
  problem->setParameter("PathOptimization/TOPPRA/interpolationMethod", hc::Parameter(interpolationMethod));
  problem->setParameter("PathOptimization/TOPPRA/gridpointMethod", hc::Parameter(gridpointMethod));
  problem->setParameter("PathOptimization/TOPPRA/N", hc::Parameter(hp::size_type(500)));
  problemParameterToCatch2Info(problem);

  CAPTURE(accLimits, *inputPath);

  auto opt = hc::pathOptimization::TOPPRA::create(problem);
  auto outputPath = opt->optimize(inputPath);

  check_path_velocity(inputPath);

  hp::vector_t q (device->configSize()), v(device->numberDof()), a(device->numberDof());
  hp::value_type
    t0 = outputPath->timeRange().first,
    t1 = outputPath->timeRange().second;

  // Check that initial and final velocities are zero
  for (auto t : { t0, t1 }) {
    outputPath->derivative(v, t, 1);
    CAPTURE(t, v);
    CHECK(v.squaredNorm() < 1e-8);
  }
  int N = 1000;
  // TODO
  // Path::timeParameterization is protected...
  // see https://github.com/humanoid-path-planner/hpp-core/pull/305
  // auto param = outputPath->timeParameterization();
  auto param = opt->lastTimeParameterization_;
  for (int i = 0; i <= N; ++i) {
    hc::value_type t = t0 + hc::value_type(i) / N * (t1-t0);
    hc::value_type
      s = param->value(t),
      sd = param->derivative(t, 1),
      sdd = param->derivative(t, 2);
    inputPath->derivative(v, s, 1);
    inputPath->derivative(a, s, 2);
    CAPTURE(t, s, sd, sdd, v, a);

    outputPath->derivative(v, t, 1);
    outputPath->derivative(a, t, 2);
    CHECK(v[0] < velLimit * constraintRelativeTol);
    if (accLimits.size() > 0)
      CHECK(a[0] < accLimits[0]*constraintRelativeTol);
  }
}

TEST_CASE("Run TOPPRA optimizer") {
  hp::DevicePtr_t device = makeDevice();
  auto problem = hc::Problem::create(device);

  std::vector<hc::PathVectorPtr_t> paths {
    makeCubicSpline(problem),
    makeCubicSpline(problem),
    makeQuinticSpline(problem, 0.2, 0.2),
    makeQuinticSpline(problem, 0.0, 0.1),
  };
  std::vector<std::string> descs {
    "cubic spline without acceleration limits",
    "cubic spline with acceleration limits",
    "quintic spline with non-zero velocity",
    "quintic spline with zero velocity at start",
  };
  std::vector<hc::vector_t> accLimitValues {
    hc::vector_t(), // No acceleration limits for cubic spline
    (hc::vector_t(1) << accLimit).finished(),
    (hc::vector_t(1) << accLimit).finished(),
    (hc::vector_t(1) << accLimit).finished(),
  };

  auto iPath = GENERATE(0, 1, 2);
  INFO(descs[iPath]);

  hp::vector_t accLimits = accLimitValues[iPath];
  CAPTURE(iPath, accLimits);
  problem->setParameter("PathOptimization/TOPPRA/accelerationLimits", hc::Parameter(accLimits));
  problem->setParameter("PathOptimization/TOPPRA/interpolationMethod", hc::Parameter(std::string("constant_acceleration")));
  problem->setParameter("PathOptimization/TOPPRA/N", hc::Parameter(hp::size_type(500)));

  auto inputPath = paths[iPath];
  auto opt = hc::pathOptimization::TOPPRA::create(problem);
  auto outputPath = opt->optimize(inputPath);

  check_path_velocity(inputPath);

  hp::vector_t q (device->configSize()), v(device->numberDof()), a(device->numberDof());
  hp::value_type
    t0 = outputPath->timeRange().first,
    t1 = outputPath->timeRange().second;

  // Check that initial and final velocities are zero
  outputPath->derivative(v, t0, 1);
  CHECK(v.isZero());
  outputPath->derivative(v, t1, 1);
  CHECK(v.isZero());
  int N = 1000;
  // TODO
  // Path::timeParameterization is protected...
  // see https://github.com/humanoid-path-planner/hpp-core/pull/305
  // auto param = outputPath->timeParameterization();
  auto param = opt->lastTimeParameterization_;
  for (int i = 0; i <= N; ++i) {
    hc::value_type t = t0 + hc::value_type(i) / N * (t1-t0);
    hc::value_type
      s = param->value(t),
      sd = param->derivative(t, 1),
      sdd = param->derivative(t, 2);
    inputPath->derivative(v, s, 1);
    inputPath->derivative(a, s, 2);
    CAPTURE(t, s, sd, sdd, v, a);

    outputPath->derivative(v, t, 1);
    outputPath->derivative(a, t, 2);
    CHECK(v[0] < velLimit * constraintRelativeTol);
    if (accLimits.size() > 0)
      CHECK(a[0] < accLimits[0]*constraintRelativeTol);
  }
}
