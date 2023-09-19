#include"tests.hh"

#include<hpp/core/straight-path.hh>

TEST_CASE("Path of length 0") {
  hp::DevicePtr_t device = makeDevice();
  auto problem = hc::Problem::create(device);

  auto addPathBefore = GENERATE(false, true);
  auto addPathAfter = GENERATE(false, true);

  CAPTURE(addPathBefore, addPathAfter);

  hc::PathVectorPtr_t pv = hc::PathVector::create(device->configSize(), device->numberDof());
  hp::vector_t q1 (device->configSize()),
               q2 (device->configSize()),
               q3 (device->configSize());
  q1 << 0.;
  q2 << 1.;
  q3 << 2.;
  if (addPathBefore)
    pv->appendPath(hc::StraightPath::create(device, q1, q2, 1.));
  pv->appendPath(hc::StraightPath::create(device, q2, q2, 0.));
  if (addPathAfter)
    pv->appendPath(hc::StraightPath::create(device, q2, q3, 1.));

  auto opt = ht::pathOptimization::TOPPRA::create(problem);
  hc::PathVectorPtr_t res;
  CHECK_NOTHROW(res = opt->optimize(pv));
  if (!addPathBefore && !addPathAfter)
    CHECK(res->length() == 0.);
  INFO("length: " << res->length());
}
