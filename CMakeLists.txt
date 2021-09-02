cmake_minimum_required(VERSION 3.3)  # CMake version check
project(FairSyn VERSION 1.0)
set(CMAKE_CXX_STANDARD 11)            # Enable c++11 standard
set(CMAKE_CXX_FLAGS "-O0 -Wall -Wextra -std=c++11 -Wfatal-errors -g -ldl --coverage -pthread") #-lgcov?
set(CMAKE_VERBOSE_MAKEFILE ON) # show everything
message("CMAKELIST: Root")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0")

set(FAIRSYN_ROOT ${PROJECT_SOURCE_DIR})
message("FAIRSYN_ROOT:" + ${FAIRSYN_ROOT})
set(ENV_VARIABLES env.sh)

set(CUDD_PATH "/home/mrychlicki/Desktop/cudd")
set(MASCOT_ROOT "/home/mrychlicki/Desktop/mascot-sds")
set(CUDD_LIBRARY_PATH "/home/mrychlicki/libraries/lib")
set(LD_LIBRARY_PATH "/home/mrychlicki/libraries/lib")
set(GOOGLE_TEST_PATH "/home/mrychlicki/Desktop/googletest")
set(SYLVAN_PATH "/home/mrychlicki/libraries")

set(CMAKE_BINARY_DIR ${FAIRSYN_ROOT}/build)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

include_directories(${FAIRSYN_ROOT}/ubdd)
include_directories(
        ${SYLVAN_PATH}/lib
        ${SYLVAN_PATH}/include)
include_directories(
        ${CUDD_PATH}
        ${CUDD_PATH}/include
        ${CUDD_PATH}/cudd
        ${CUDD_PATH}/lib
        ${CUDD_PATH}/util
        ${CUDD_PATH}/dddmph)

include_directories(${CUDD_LIBRARY_PATH})
#include_directories(
#        ${HOAF_LIB_PATH}/
#        ${HOAF_LIB_PATH}/cpphoafparser
#        ${HOAF_LIB_PATH}/cpphoafparser/consumer
#        ${HOAF_LIB_PATH}/cpphoafparser/parser
#        ${HOAF_LIB_PATH}/cpphoafparser/util
#        ${HOAF_LIB_PATH}/cpphoafparser/ast)

add_subdirectory(${FAIRSYN_ROOT}/ubdd)
add_subdirectory(${FAIRSYN_ROOT}/tests)
