message(CMakeLists:FairSyn/cmake/FindCppHoafParser.cmake)
find_package(PkgConfig)
pkg_check_modules(PC_CppHoafParser QUIET CppHoafParser)

if (CHP_HOME)
    find_path(CppHoafParser_INCLUDE_DIR
            NAMES parser/hoa_parser.hh consumer/hoa_consumer.hh
            PATHS ${CHP_HOME}
            PATH_SUFFIXES cpphoafparser
            NO_DEFAULT_PATH
            )
else ()
    find_path(CppHoafParser_INCLUDE_DIR
            NAMES parser/hoa_parser.hh consumer/hoa_consumer.hh
            PATH_SUFFIXES cpphoafparser
            )
endif ()

set(CppHoafParser_VERSION ${PC_CppHoafParser_VERSION})

mark_as_advanced(CppHoafParser_FOUND CppHoafParser_INCLUDE_DIR CppHoafParser_VERSION)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CppHoafParser
        REQUIRED_VARS CppHoafParser_INCLUDE_DIR
        VERSION_VAR CppHoafParser_VERSION
        )

if (CppHoafParser_FOUND AND NOT TARGET CppHoafParser::CppHoafParser)
    add_library(CppHoafParser::CppHoafParser INTERFACE IMPORTED)
    set_target_properties(CppHoafParser::CppHoafParser PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${CppHoafParser_INCLUDE_DIR}/consumer;${CppHoafParser_INCLUDE_DIR}/parser;${CppHoafParser_INCLUDE_DIR}/ast;${CppHoafParser_INCLUDE_DIR}/util"
            )
endif ()
#message(CppHoafParser_INCLUDE_DIR: ${CppHoafParser_INCLUDE_DIR})