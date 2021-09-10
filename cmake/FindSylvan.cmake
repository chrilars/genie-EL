message(mascot/fairsyn/cmake/FindSylvan.cmake)

find_package(PkgConfig)
pkg_check_modules(PC_Sylvan QUIET Sylvan)

find_path(Sylvan_INCLUDE_DIR NAMES sylvan.h)
find_library(Sylvan_STATIC_LIB_DIR NAMES libSylvan.a)
find_library(Sylvan_SHARED_LIB_DIR NAMES libSylvan.so)

set(Sylvan_VERSION ${PC_Sylvan_VERSION})

mark_as_advanced(Sylvan_FOUND Sylvan_INCLUDE_DIR Sylvan_VERSION)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Sylvan
        REQUIRED_VARS Sylvan_INCLUDE_DIR
        VERSION_VAR Sylvan_VERSION
        )

if(Sylvan_FOUND)
    #Set include dirs to parent, to enable includes like #include <Sylvan/document.h>
    get_filename_component(Sylvan_INCLUDE_DIRS ${Sylvan_INCLUDE_DIR} DIRECTORY)
endif()

if(Sylvan_FOUND AND NOT TARGET Sylvan::Sylvan)
    add_library(Sylvan::Sylvan STATIC IMPORTED)
    set_target_properties(
            Sylvan::Sylvan
                PROPERTIES
                    IMPORTED_LOCATION ${Sylvan_STATIC_LIB_DIR}
                    INTERFACE_INCLUDE_DIRECTORIES "${Sylvan_INCLUDE_DIRS}"
            )
endif()
#add_link_options(--coverage)
