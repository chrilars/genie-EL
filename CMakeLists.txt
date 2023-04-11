#  CMakeList.txt
#
#  created on: 10.09.2021 (In Progress)
#      author: Mateusz Rychlicki
#   Based on: https://github.com/pabloariasal/modern-cmake-sample/
#

cmake_minimum_required(VERSION 3.13...3.19 FATAL_ERROR)
project(libgenie VERSION 1.0.0 LANGUAGES CXX)
message("CMakeLists:Genie")

list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

##############################################
# Find system dependencies
find_package(cudd)

find_package(sylvan)

find_package(CppHoafParser)
##############################################
# Create target and set properties

add_library(genie
        src/ubdd/CuddUBDD.cpp
        src/ubdd/CuddUBDDMintermIterator.cpp
        src/ubdd/SylvanUBDD.cpp
        src/ubdd/SylvanUBDDMintermIterator.cpp
        include/lib/BaseRabinAutomaton.hh
        include/lib/BaseFixpoint.hh
        include/lib/ParallelRabinRecurse.hh)


#Add an alias so that library can be used inside the build tree, e.g. when testing
add_library(Genie::genie ALIAS genie)

#Set target properties
target_include_directories(genie
        PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        )

target_compile_features(genie PRIVATE cxx_auto_type)

target_link_libraries(genie ${cudd_LIBRARIES} ${sylvan_LIBRARIES} CppHoafParser)
target_include_directories(genie PUBLIC ${cudd_INCLUDE_DIRS})
target_compile_options(genie PUBLIC ${cudd_CFLAGS_OTHER})

if(NOT DEFINED CMAKE_BUILD_TYPE)
    message(STATUS "No mode selected. Using release mode as default.")
    set(CMAKE_BUILD_TYPE "release")
endif()

set(CMAKE_VERBOSE_MAKEFILE ON) # show everything
add_link_options(-std=c++11 -ldl -O3)
if(CMAKE_BUILD_TYPE MATCHES debug)
    message(STATUS "Debug mode.")
    add_link_options(-Wfatal-errors -Wall -Wextra)
    target_compile_options(genie PRIVATE
            $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:
            -Wall -Wextra>)
elseif(CMAKE_BUILD_TYPE MATCHES release)
    message(STATUS "Release mode.")
    target_compile_options(genie PRIVATE
            $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>: >)
else()
    message(STATUS "Undefined mode.")
    target_compile_options(genie PRIVATE
            $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>: >)
endif()

##############################################
# Installation instructions

include(GNUInstallDirs)
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/Genie)

install(TARGETS genie
        EXPORT genie-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )

#This is required so that the exported target has the name Genie and not genie
set_target_properties(genie PROPERTIES EXPORT_NAME Genie)

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

#Export the targets to a script
install(EXPORT genie-targets
        FILE
        GenieTargets.cmake
        NAMESPACE
        Genie::
        DESTINATION
        ${INSTALL_CONFIGDIR}
        )

#Create a ConfigVersion.cmake file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/GenieConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion
)

configure_package_config_file(${CMAKE_CURRENT_LIST_DIR}/cmake/GenieConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/GenieConfig.cmake
        INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
        )

#Install the config, configversion and custom find modules
install(FILES
        ${CMAKE_CURRENT_LIST_DIR}/cmake/Findcudd.cmake
        ${CMAKE_CURRENT_LIST_DIR}/cmake/Findsylvan.cmake
        ${CMAKE_CURRENT_LIST_DIR}/cmake/FindCppHoafParser.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/GenieConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/GenieConfigVersion.cmake
        DESTINATION ${INSTALL_CONFIGDIR}
        )

##############################################
## Exporting from the build tree
configure_file(${CMAKE_CURRENT_LIST_DIR}/cmake/Findcudd.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/Findcudd.cmake
        COPYONLY)

configure_file(${CMAKE_CURRENT_LIST_DIR}/cmake/Findsylvan.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/Findsylvan.cmake
        COPYONLY)

configure_file(${CMAKE_CURRENT_LIST_DIR}/cmake/FindCppHoafParser.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/FindCppHoafParser.cmake
        COPYONLY)

export(EXPORT genie-targets
        FILE ${CMAKE_CURRENT_BINARY_DIR}/GenieTargets.cmake
        NAMESPACE Genie::)

#Register package in the User Package Registry
export(PACKAGE Genie)
