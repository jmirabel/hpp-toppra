find_package(pybind11 REQUIRED)

message("Found ${PYTHON_INCLUDE_DIR} ${PYTHON_EXECUTABLE}")

pybind11_add_module(hpp_toppra_cpp
  path.cc
)
target_link_libraries(hpp_toppra_cpp
  PUBLIC
  hpp-core::hpp-core
  hpp-manipulation::hpp-manipulation
  toppra::toppra
)
