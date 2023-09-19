#include"catch.hpp"

#include<hpp/toppra/toppra.hh>

#include<hpp/core/problem.hh>
#include<hpp/core/path-vector.hh>
#include<hpp/core/path/spline.hh>
#include<hpp/core/steering-method/spline.hh>

namespace hc = hpp::core;
namespace hp = hpp::pinocchio;
namespace ht = hpp::toppra;

typedef hc::steeringMethod::Spline< hc::path::BernsteinBasis, 1 > SteerSplineOrder1;
typedef hc::steeringMethod::Spline< hc::path::BernsteinBasis, 3 > SteerSplineOrder3;
typedef hc::steeringMethod::Spline< hc::path::BernsteinBasis, 5 > SteerSplineOrder5;

hp::DevicePtr_t makeDevice();

hc::PathVectorPtr_t makeCubicSpline(hc::ProblemPtr_t p);
