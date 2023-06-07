cmake_minimum_required(VERSION 3.4)

if (NOT $ENV{CONDA_PREFIX} STREQUAL "")
    message(STATUS "CURRENT CONDA PREFIX: $ENV{CONDA_PREFIX}")
    include_directories($ENV{CONDA_PREFIX}/include)
    link_directories($ENV{CONDA_PREFIX}/lib)
endif ()

set(CMAKE_CXX_COMPILE g++)
set(CMAKE_C_COMPILE gcc)

project(ExtractHAIRs)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE DEBUG)
set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXX_FLAGS} -O0 -Wall -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 ")

set(CMAKE_INCLUDE_CURRENT_DIR ON)
#add_subdirectory(lib/bamtools-2.3.0)
#include_directories (lib/bamtools-2.3.0/src)
#include_directories(lib/tclap-1.2.1/include)
aux_source_directory(. EXTRACT_SRC)
aux_source_directory(./htslib HTS_SRC)
aux_source_directory(./samtools SAM_SRC)

#set(LIBSAM libbam.a)
set(LIBHTS libhts.a)


add_executable(ExtractHAIRs ${EXTRACT_SRC})

target_link_libraries(ExtractHAIRs -lhts z)
