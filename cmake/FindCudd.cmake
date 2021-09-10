message(mascot/fairsyn/cmake/FindCudd.cmake)
find_package(PkgConfig)
pkg_check_modules(PC_Cudd QUIET Cudd)

find_path(Cudd_INCLUDE_DIR NAMES cudd.h uddObj.hh dddmp.h)
find_library(Cudd_STATIC_LIB_DIR NAMES libcudd.a)
find_library(Cudd_SHARED_LIB_DIR NAMES libcudd.so)

set(Cudd_VERSION ${PC_Cudd_VERSION})

mark_as_advanced(Cudd_FOUND Cudd_INCLUDE_DIR Cudd_VERSION)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Cudd
        REQUIRED_VARS Cudd_INCLUDE_DIR
        VERSION_VAR Cudd_VERSION
        )

if(Cudd_FOUND)
    #Set include dirs to parent, to enable includes like #include <rapidjson/document.h>
    get_filename_component(Cudd_INCLUDE_DIRS ${Cudd_INCLUDE_DIR} DIRECTORY)
endif()

if(Cudd_FOUND AND NOT TARGET Cudd::Cudd)
    add_library(Cudd::Cudd STATIC IMPORTED)
    set_target_properties(
            Cudd::Cudd
                PROPERTIES
                    IMPORTED_LOCATION ${Cudd_STATIC_LIB_DIR}
                    INTERFACE_INCLUDE_DIRECTORIES "${Cudd_INCLUDE_DIRS}"
            )
endif()