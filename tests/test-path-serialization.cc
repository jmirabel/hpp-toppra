#include <hpp/core/straight-path.hh>
#include <hpp/manipulation/device.hh>
#include <hpp/pinocchio/serialization.hh>
#include <hpp/pinocchio/urdf/util.hh>

#include "../src/serialization.hh"
#include "tests.hh"

TEST_CASE("serialize and read") {
  std::string pathStringVersion;
  {
    hp::DevicePtr_t device = makeDevice();
    auto problem = hc::Problem::create(device);

    auto path = makeCubicSpline(problem);
    std::ostringstream oss;
    oss << *path;
    pathStringVersion = oss.str();

    std::ofstream fs("path.txt");
    hpp::serialization::text_oarchive archive(fs);
    archive.insert(device->name(), device.get());
    archive << path;
  }
  INFO(pathStringVersion);
  {
    hp::DevicePtr_t device = makeDevice();

    std::ifstream fs("path.txt");
    hpp::serialization::text_iarchive archive(fs);
    archive.insert(device->name(), device.get());
    hc::PathVectorPtr_t path;
    archive >> path;

    std::ostringstream oss;
    oss << *path;
    std::string pathAfterSerialization = oss.str();
    INFO(pathAfterSerialization);
    CHECK(pathStringVersion == pathAfterSerialization);
  }
}

/*
TEST_CASE("inspect serialized path") {
  hpp::manipulation::Device::create("fake_world");
  hp::DevicePtr_t device = hp::Device::create("world");
  hpp::pinocchio::urdf::loadModel(
      device, 0, "", "anchor", "package://ur_description_eureka/urdf/ur3e.urdf",
      "package://ur_description_eureka/srdf/ur3e.srdf");
  hc::PathPtr_t path;

  hc::parser::serializePath<hpp::serialization::binary_iarchive>(
      device, path, "./last-hpp-toppra-path.bin");
  std::cout << *path << std::endl;
}
*/
