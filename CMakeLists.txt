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
set(PROJECT_EXPORT_NO_TARGET TRUE)

set(PROJECT_NAME hpp-toppra)
set(PROJECT_DESCRIPTION "Bridge between HPP and TOPPRA.")

include(cmake/hpp.cmake)
include(cmake/boost.cmake)

compute_project_args(PROJECT_ARGS LANGUAGES CXX)
project(${PROJECT_NAME} ${PROJECT_ARGS})

set(CMAKE_CXX_STANDARD 14)
add_project_dependency(hpp-core REQUIRED)
add_project_dependency(hpp-manipulation REQUIRED)
add_project_dependency(toppra REQUIRED)

include(${HPP_CORE_CMAKE_PLUGIN})

message(STATUS "Found toppra ${toppra_VERSION}")
hpp_add_plugin(toppra
  EXPORT toppraTargets
  SOURCES src/toppra
  LINK_DEPENDENCIES PUBLIC hpp-core::hpp-core toppra::toppra)

add_subdirectory(python)

add_executable(unittests
  tests/tests.cc
  tests/test-ill-formed-paths.cc
  tests/test-path-serialization.cc
  src/toppra.cc
)
add_test(unittests unittests)
target_link_libraries(unittests hpp-core::hpp-core hpp-manipulation::hpp-manipulation toppra::toppra)# Catch2::Catch2)

project_install_component(toppra)
install(FILES src/toppra.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/hpp/core/path-optimization/)
