message(CMakeLists:FairSyn/cmake/Findcudd.cmake)

if (cudd_HOME)
    message("cudd home")#debug
    find_path(cudd_INCLUDE_DIR cudd.h PATHS "${cudd_HOME}/include" NO_DEFAULT_PATH)
    find_library(cudd_LIBRARY libcudd.a cudd PATHS "${cudd_HOME}/lib" NO_DEFAULT_PATH)
else ()
    message("cudd no home")#debug
    find_path(cudd_INCLUDE_DIR cudd.h)
    find_library(cudd_LIBRARY libcudd.a cudd)
endif ()

message("${cudd_LIBRARY}")#debug
message("${cudd_INCLUDE_DIR}")#debug

set(cudd_LIBRARIES ${cudd_LIBRARY})
set(cudd_INCLUDE_DIRS ${cudd_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(cudd DEFAULT_MSG cudd_LIBRARY cudd_INCLUDE_DIR)

mark_as_advanced(cudd_INCLUDE_DIR cudd_LIBRARY)
