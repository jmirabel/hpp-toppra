#
# Copyright (c) 2022 Eureka Robotics
# Authors: Joseph Mirabel
#
#

cmake_minimum_required(VERSION 3.1)

set(CXX_DISABLE_WERROR TRUE)
set(DOXYGEN_USE_MATHJAX YES)
set(DOXYGEN_USE_TEMPLATE_CSS TRUE)
set(PROJECT_USE_CMAKE_EXPORT TRUE)

set(PROJECT_NAME hpp-toppra)
set(PROJECT_DESCRIPTION "Bridge between HPP and TOPPRA.")

include(cmake/hpp.cmake)
include(cmake/boost.cmake)

compute_project_args(PROJECT_ARGS LANGUAGES CXX)
project(${PROJECT_NAME} ${PROJECT_ARGS})

set(CMAKE_CXX_STANDARD 14)
add_project_dependency(hpp-core REQUIRED)
add_project_dependency(hpp-manipulation REQUIRED)
add_project_dependency(toppra 0.6.2 REQUIRED)

include(${HPP_CORE_CMAKE_PLUGIN})

set(${PROJECT_NAME}_HEADERS
  include/hpp/toppra/toppra.hh
)

add_library(${PROJECT_NAME}
  SHARED
  src/toppra.cc
  ${${PROJECT_NAME}_HEADERS}
)
target_link_libraries(hpp-toppra PUBLIC hpp-core::hpp-core toppra::toppra)
target_include_directories(hpp-toppra PUBLIC $<INSTALL_INTERFACE:include>)
install(
  TARGETS ${PROJECT_NAME}
  EXPORT ${TARGETS_EXPORT_NAME}
  DESTINATION lib)

hpp_add_plugin(toppra
  EXPORT toppraTargets
  SOURCES
  src/piecewise-polynomial.hh
  src/serialization.hh
  src/toppra-plugin.cc
  LINK_DEPENDENCIES PUBLIC hpp-toppra)
#target_link_libraries(toppra PUBLIC hpp-toppra)

message(STATUS "Found toppra ${toppra_VERSION}")

add_subdirectory(python)

add_executable(unittests
  tests/tests.cc
  tests/test-ill-formed-paths.cc
  tests/test-path-serialization.cc
)
add_test(unittests unittests)
target_link_libraries(unittests hpp-manipulation::hpp-manipulation toppra::toppra hpp-toppra)# Catch2::Catch2)

project_install_component(toppra)
