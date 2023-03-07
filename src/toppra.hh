#ifndef HPP_CORE_PATH_OPTIMIZATION_TOPPRA_HH
#define HPP_CORE_PATH_OPTIMIZATION_TOPPRA_HH

#include <hpp/core/path-optimizer.hh>

#include <toppra/toppra.hpp>

namespace hpp {
namespace core {
namespace pathOptimization {

class TOPPRA;
typedef shared_ptr<TOPPRA> TOPPRAPtr_t;

class TOPPRA : public PathOptimizer {
 public:
   enum InterpolationMethod {
     ConstantAcceleration,
     Hermite,
   };
   enum GridpointMethod {
     EvenlyTimeSpaced,
     EvenlyParamSpaced,
   };

  static TOPPRAPtr_t create(const ProblemConstPtr_t &p) {
    return TOPPRAPtr_t(new TOPPRA(p));
  }

  PathVectorPtr_t optimize(const PathVectorPtr_t &path);

  // TODO remove when
  // https://github.com/humanoid-path-planner/hpp-core/pull/305
  // is released.
  TimeParameterizationPtr_t lastTimeParameterization_;

 protected:
  using PathOptimizer::PathOptimizer;

 private:
  void inputSerialization(PathPtr_t path) const;
  toppra::LinearConstraintPtrs constraints();
  InterpolationMethod interpolationMethod() const;
  GridpointMethod gridpointMethod() const;
};  // class TOPPRA

}  // namespace pathOptimization
}  // namespace core
}  // namespace hpp

#endif  // HPP_CORE_PATH_OPTIMIZATION_TOPPRA_HH
