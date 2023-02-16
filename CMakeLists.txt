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

add_project_dependency(hpp-core REQUIRED)
add_project_dependency(toppra REQUIRED)

include(${HPP_CORE_CMAKE_PLUGIN})

message(STATUS "Found toppra ${toppra_VERSION}")
hpp_add_plugin(toppra
  EXPORT toppraTargets
  SOURCES src/toppra
  LINK_DEPENDENCIES PUBLIC hpp-core::hpp-core toppra::toppra)

if(${CMAKE_VERSION} VERSION_GREATER 3.11.4)
  Include(FetchContent)

  FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v3.0.1)

  FetchContent_MakeAvailable(Catch2)
else()
  include(ExternalProject)
  find_package(Git REQUIRED)

  ExternalProject_Add(
          catch2
          PREFIX ${CMAKE_BINARY_DIR}/catch2
          GIT_REPOSITORY https://github.com/catchorg/Catch2.git
          GIT_TAG v3.0.1
          GIT_SHALLOW TRUE
          TIMEOUT 10
          UPDATE_COMMAND ${GIT_EXECUTABLE} pull
          CONFIGURE_COMMAND ""
          BUILD_COMMAND ""
          INSTALL_COMMAND ""
          LOG_DOWNLOAD ON
          )

  ExternalProject_Get_Property(catch2 SOURCE_DIR)
  set(CATCH_INCLUDE_DIR ${SOURCE_DIR}/single_include)
  message("CATCH_INCLUDE_DIR: ${CATCH_INCLUDE_DIR}")
  add_library(Catch2 INTERFACE)
  target_include_directories(Catch2 INTERFACE ${CATCH_INCLUDE_DIR})

  add_library(Catch2::Catch2 ALIAS Catch2)
endif()


add_executable(unittests
  tests/tests.cc
  src/toppra.cc
)
add_test(unittests unittests)
target_link_libraries(unittests hpp-core::hpp-core toppra::toppra Catch2::Catch2WithMain)

project_install_component(toppra)
install(FILES src/toppra.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/hpp/core/path-optimization/)
