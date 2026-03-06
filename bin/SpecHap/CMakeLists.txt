cmake_minimum_required(VERSION 3.6)

project(SpecHap)

set(CMAKE_CXX_STANDARD 11)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

if(CMAKE_PREFIX_PATH)
  link_directories(${CMAKE_PREFIX_PATH}/lib)
endif()

find_package(htslib REQUIRED)
include_directories(${HTSlib_INCLUDE_DIRS})

# arpack is a Fortran library with no headers - find the library directly
find_library(ARPACK_LIB NAMES arpack PATHS ${CMAKE_PREFIX_PATH}/lib $ENV{CONDA_PREFIX}/lib)
if(NOT ARPACK_LIB)
  message(FATAL_ERROR "arpack library not found")
endif()
message(STATUS "   ARPACK library: ${ARPACK_LIB}")

aux_source_directory(. SPEC_SRC)

add_executable(SpecHap ${SPEC_SRC})
target_link_libraries(SpecHap ${HTSlib_LIBRARIES} ${ARPACK_LIB})

install(TARGETS SpecHap DESTINATION bin)
