message(CMakeLists:FairSyn/cmake/Findcudd.cmake)

find_path(cudd_INCLUDE_DIR cudd.h)
find_library(cudd_LIBRARY libcudd.a cudd)

set(cudd_LIBRARIES ${cudd_LIBRARY})
set(cudd_INCLUDE_DIRS ${cudd_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(cudd DEFAULT_MSG cudd_LIBRARY cudd_INCLUDE_DIR)

mark_as_advanced(cudd_INCLUDE_DIR cudd_LIBRARY)
