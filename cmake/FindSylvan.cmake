message(CMakeLists:FairSyn/cmake/FindSylvan.cmake)
find_package(PkgConfig)
pkg_check_modules(PC_Sylvan QUIET Sylvan)

find_path(Sylvan_INCLUDE_DIR sylvan.h)
find_library(Sylvan_LIBRARY libSylvan.a)
if (NOT Sylvan_LIBRARY)
    # Name of this project has been changed from "Sylvan" to "sylvan"
    find_library(Sylvan_LIBRARY libsylvan.a)
endif ()

set(Sylvan_LIBRARIES ${Sylvan_LIBRARY})
set(Sylvan_INCLUDE_DIRS ${Sylvan_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Sylvan DEFAULT_MSG Sylvan_LIBRARY Sylvan_INCLUDE_DIR)

mark_as_advanced(Sylvan_INCLUDE_DIR Sylvan_LIBRARY)