#include"catch.hpp"

#include<hpp/core/problem.hh>
#include<hpp/core/path-vector.hh>
#include<hpp/core/path/spline.hh>
#include<hpp/core/steering-method/spline.hh>

#include"../src/toppra.hh"

namespace hc = hpp::core;
namespace hp = hpp::pinocchio;

typedef hc::steeringMethod::Spline< hc::path::BernsteinBasis, 1 > SteerSplineOrder1;
typedef hc::steeringMethod::Spline< hc::path::BernsteinBasis, 3 > SteerSplineOrder3;
typedef hc::steeringMethod::Spline< hc::path::BernsteinBasis, 5 > SteerSplineOrder5;

hp::DevicePtr_t makeDevice();

hc::PathVectorPtr_t makeCubicSpline(hc::ProblemPtr_t p);
