cmake_minimum_required (VERSION 2.6)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")

project (cylindrical_fdtd)

find_package(HDF5 REQUIRED COMPONENTS C CXX)
include_directories(${HDF5_INCLUDE_DIRS})

add_executable(cylindrical_fdtd cylindrical_fdtd.cpp)
target_link_libraries(cylindrical_fdtd ${HDF5_LIBRARIES})
