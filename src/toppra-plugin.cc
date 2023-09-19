#include <hpp/toppra/toppra.hh>

#include <hpp/core/plugin.hh>
#include <hpp/core/problem-solver.hh>

namespace hpp {
namespace toppra {

using hpp::core::ProblemSolverPlugin;
using hpp::core::ProblemSolverPtr_t;

class TOPPRAPlugin : public ProblemSolverPlugin {
 public:
  TOPPRAPlugin() : ProblemSolverPlugin("TOPPRAPlugin", "0.0") {}

 protected:
  virtual bool impl_initialize(ProblemSolverPtr_t ps) {
    ps->pathOptimizers.add("TOPPRA", pathOptimization::TOPPRA::create);
    return true;
  }
};
}  // namespace toppra
}  // namespace hpp

HPP_CORE_DEFINE_PLUGIN(hpp::toppra::TOPPRAPlugin)
