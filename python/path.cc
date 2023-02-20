#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

#include<sstream>
#include<hpp/pinocchio/device.hh>
#include<hpp/pinocchio/urdf/util.hh>
#include<hpp/core/path.hh>
#include<hpp/core/path-vector.hh>
#include<hpp/manipulation/device.hh>

#include <toppra/geometric_path.hpp>

#include "../src/serialization.hh"

namespace py = pybind11;

namespace hpp {
namespace pinocchio {
void expose(py::module &m) {
  py::class_<Device, DevicePtr_t>(m, "Device")
    .def_static("create", &Device::create)
    ;

  {
    py::module sm = m.def_submodule("urdf");
    sm.def("loadModel", [](
          const DevicePtr_t& robot, const FrameIndex& baseFrame,
          const std::string& prefix, const std::string& rootType,
          const std::string& urdfPath, const std::string& srdfPath
          // , const SE3& bMr = SE3::Identity()
          ) {
        urdf::loadModel(robot, baseFrame, prefix, rootType, urdfPath, srdfPath, SE3::Identity());
      });
  }
}
}
namespace core {

class TOPPRAPathWrapper : public toppra::GeometricPath {
 public:
  TOPPRAPathWrapper(PathPtr_t path)
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

PathPtr_t readBinPath(DevicePtr_t device, std::string filename) {
  PathPtr_t path;
  parser::serializePath<serialization::binary_iarchive>(device, path, filename);
  return path;
}

void expose(py::module &m) {
  py::class_<Path, PathPtr_t>(m, "Path")
    .def("eval", [](const PathPtr_t& path, value_type time) {
        bool success;
        Configuration_t q = path->eval(time, success);
        return py::make_tuple(q, success);
    })
    .def("derivative", [](const PathPtr_t& path, value_type time, size_type order) {
        vector_t der (path->outputDerivativeSize());
        path->derivative(der, time, order);
        return der;
    })
    .def("velocityBound", [](const PathPtr_t& path, value_type t0, value_type t1) {
        vector_t res (path->outputDerivativeSize());
        path->velocityBound(res, t0, t1);
        return res;
    })
    .def_property_readonly("timeRange", [](Path const& p) { return p.timeRange(); })
    .def("length", &Path::length)
    .def("str", [](const Path& p) {
        std::ostringstream oss;
        oss << p;
        return oss.str();
      })
    ;
  py::class_<PathVector, PathVectorPtr_t, Path>(m, "PathVector")
    .def("numberPaths", &PathVector::numberPaths)
    .def("pathAtRank", &PathVector::pathAtRank)
  ;
  py::class_<TOPPRAPathWrapper,
    std::shared_ptr<TOPPRAPathWrapper>,
    toppra::GeometricPath>(m, "TOPPRAPathWrapper")
    .def(py::init<PathPtr_t>())
  ;

  m.def("readBinPath", &readBinPath);
}
}
}

PYBIND11_MODULE(hpp_toppra_cpp, m)
{
  py::module_::import("toppra.cpp");
  {
    py::module sm = m.def_submodule("pinocchio");
    hpp::pinocchio::expose(sm);
  }
  {
    py::module sm = m.def_submodule("core");
    hpp::core::expose(sm);
  }

  hpp::manipulation::Device::create("fake");
}
