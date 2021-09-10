#
# CuddUBDD.hh
#
#  created on: 10.09.2021 (In Progress)
#      author: Mateusz Rychlicki
#   Based on: https://github.com/pabloariasal/modern-cmake-sample/
#

cmake_minimum_required(VERSION 3.13...3.19 FATAL_ERROR)
project(libfairsyn VERSION 1.0.0 LANGUAGES CXX)
message(mascot/farisyn)

#Make sure that custom modules like FindRapidJSON are found
list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

##############################################
# Find system dependencies
set(MIN_Cudd_VERSION 3.0)
find_package(Cudd ${MIN_Cudd_VERSION} REQUIRED)

set(MIN_Sylvan_VERSION 1.6)
find_package(Sylvan REQUIRED)

set(MIN_CppHoafParser_VERSION 2.0)
find_package(CppHoafParser REQUIRED ${MIN_CppHoafParser_VERSION})
##############################################
# Create target and set properties

add_library(fairsyn
        src/ubdd/CuddUBDD.cpp
        src/ubdd/CuddUBDDMintermIterator.cpp
        src/ubdd/SylvanUBDD.cpp
        src/ubdd/SylvanUBDDMintermIterator.cpp
#        src/lib/Arena.cpp
#        src/lib/FixedPoint.cpp
        src/utils/hoa_consumer_build_rabin.cpp
        src/utils/TicToc.cpp
        )


#Add an alias so that library can be used inside the build tree, e.g. when testing
add_library(FairSyn::fairsyn ALIAS fairsyn)

#Set target properties
target_include_directories(fairsyn
        PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        )

target_compile_features(fairsyn PRIVATE cxx_auto_type)
target_compile_options(fairsyn PRIVATE
        $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:
        -Wall -Wextra -Wpedantic>)

target_link_libraries(fairsyn PRIVATE Cudd::Cudd Sylvan::Sylvan CppHoafParser::CppHoafParser)
add_link_options(-std=c++11 -g -ldl) # Mascot
add_link_options(-O0 -Wall -Wextra -Wfatal-errors) # Debug
add_link_options(--coverage -pthread)

##############################################
# Installation instructions

include(GNUInstallDirs)
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/FairSyn)

install(TARGETS fairsyn
        EXPORT fairsyn-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )

#This is required so that the exported target has the name FairSyn and not fairsyn
set_target_properties(fairsyn PROPERTIES EXPORT_NAME FairSyn)

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

#Export the targets to a script
install(EXPORT fairsyn-targets
        FILE
        FairSynTargets.cmake
        NAMESPACE
        FairSyn::
        DESTINATION
        ${INSTALL_CONFIGDIR}
        )

#Create a ConfigVersion.cmake file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/FairSynConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion
)

configure_package_config_file(${CMAKE_CURRENT_LIST_DIR}/cmake/FairSynConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/FairSynConfig.cmake
        INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
        )

#Install the config, configversion and custom find modules
install(FILES
        ${CMAKE_CURRENT_LIST_DIR}/cmake/FindCudd.cmake
        ${CMAKE_CURRENT_LIST_DIR}/cmake/FindSylvan.cmake
        ${CMAKE_CURRENT_LIST_DIR}/cmake/FindCppHoafParser.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/FairSynConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/FairSynConfigVersion.cmake
        DESTINATION ${INSTALL_CONFIGDIR}
        )

##############################################
## Exporting from the build tree
configure_file(${CMAKE_CURRENT_LIST_DIR}/cmake/FindCudd.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/FindCudd.cmake
        COPYONLY)

configure_file(${CMAKE_CURRENT_LIST_DIR}/cmake/FindSylvan.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/FindSylvan.cmake
        COPYONLY)

configure_file(${CMAKE_CURRENT_LIST_DIR}/cmake/FindCppHoafParser.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/FindCppHoafParser.cmake
        COPYONLY)

export(EXPORT fairsyn-targets
        FILE ${CMAKE_CURRENT_BINARY_DIR}/FairSynTargets.cmake
        NAMESPACE FairSyn::)

#Register package in the User Package Registry
export(PACKAGE FairSyn)

##############################################
## Add test
enable_testing()
add_subdirectory(test)
