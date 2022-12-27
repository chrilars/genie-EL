message(CMakeLists:FairSyn/cmake/Findsylvan.cmake)
find_package(PkgConfig)
pkg_check_modules(PC_sylvan QUIET sylvan)

# Name of this project has been changed from "sylvan" to "sylvan"
find_path(sylvan_INCLUDE_DIR sylvan.h)
find_library(sylvan_LIBRARY libsylvan.so)


set(sylvan_LIBRARIES ${sylvan_LIBRARY})
set(sylvan_INCLUDE_DIRS ${sylvan_INCLUDE_DIR})


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(sylvan DEFAULT_MSG sylvan_LIBRARY sylvan_INCLUDE_DIR)

mark_as_advanced(sylvan_INCLUDE_DIR sylvan_LIBRARY)