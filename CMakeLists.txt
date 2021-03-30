cmake_minimum_required (VERSION 3.11.4)
project (limbo C CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17) 
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
endif()

include(CheckCXXCompilerFlag)


SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

add_compile_options("-Wall")

add_subdirectory(keccak)

set(BANQUET_SRCS
  BristolCircuit.cpp
  limbo.cpp
  limbo_instances.cpp
  field.cpp
  tree.cpp
  tape.cpp
  randomness.c
  )


set(PCLMUL_FILES field.cpp)
CHECK_CXX_COMPILER_FLAG(-mpclmul COMPILER_HAS_M_PCLMUL)
CHECK_CXX_COMPILER_FLAG(-msse2 COMPILER_HAS_M_SSE2)
CHECK_CXX_COMPILER_FLAG(-msse4 COMPILER_HAS_M_SSE4)
if (COMPILER_HAS_M_PCLMUL AND COMPILER_HAS_M_SSE2 AND COMPILER_HAS_M_SSE4)
  add_compile_options(-mpclmul -msse2 -msse4)
else()
  message(
    ERROR
    "compiler does not have at least one of flag (pclmul, sse2, sse4) which are needed"  )
endif()

add_library(limbo_static STATIC ${BANQUET_SRCS})
target_link_libraries(limbo_static PUBLIC keccak) 

add_executable(bench bin/bench.cpp)
target_link_libraries(bench limbo_static)
