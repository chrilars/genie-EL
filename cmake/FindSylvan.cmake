message(CMakeLists:FairSyn/cmake/FindSylvan.cmake)
find_package(PkgConfig)
pkg_check_modules(PC_Sylvan QUIET Sylvan)

if (Sylvan_HOME)
    message("Sylvan home")#debug
    find_path(Sylvan_INCLUDE_DIR sylvan.h PATHS "${Sylvan_HOME}/include" NO_DEFAULT_PATH)
    find_library(Sylvan_LIBRARY libSylvan.a PATHS "${Sylvan_HOME}/lib" NO_DEFAULT_PATH)
else ()
    message("Sylvan no home")#debug
    find_path(Sylvan_INCLUDE_DIR sylvan.h)
    find_library(Sylvan_LIBRARY libSylvan.a)
endif ()

message("${Sylvan_LIBRARY}")#debug
message("${Sylvan_INCLUDE_DIR}")#debug

set(Sylvan_LIBRARIES ${Sylvan_LIBRARY})
set(Sylvan_INCLUDE_DIRS ${Sylvan_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Sylvan DEFAULT_MSG Sylvan_LIBRARY Sylvan_INCLUDE_DIR)

mark_as_advanced(Sylvan_INCLUDE_DIR Sylvan_LIBRARY)