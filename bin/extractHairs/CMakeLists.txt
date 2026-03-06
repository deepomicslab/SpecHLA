cmake_minimum_required(VERSION 3.5)

project(ExtractHAIRs)

set(CMAKE_CXX_STANDARD 14)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

if(CMAKE_PREFIX_PATH)
  include_directories(${CMAKE_PREFIX_PATH}/include)
  link_directories(${CMAKE_PREFIX_PATH}/lib)
endif()

set(CMAKE_INCLUDE_CURRENT_DIR ON)
aux_source_directory(. EXTRACT_SRC)
aux_source_directory(./htslib HTS_SRC)
aux_source_directory(./samtools SAM_SRC)

add_executable(ExtractHAIRs ${EXTRACT_SRC})
target_link_libraries(ExtractHAIRs -lhts z)

install(TARGETS ExtractHAIRs DESTINATION bin)
