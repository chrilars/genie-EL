message(CMakeLists:Genie/cmake/FindCppHoafParser.cmake)
pkg_check_modules(PC_CppHoafParser QUIET CppHoafParser)
find_path(CppHoafParser_INCLUDE_DIR
        NAMES parser/hoa_parser.hh consumer/hoa_consumer.hh
        PATH_SUFFIXES cpphoafparser
        )

set(CppHoafParser_VERSION ${PC_CppHoafParser_VERSION})

mark_as_advanced(CppHoafParser_FOUND CppHoafParser_INCLUDE_DIR CppHoafParser_VERSION)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CppHoafParser CppHoafParser_INCLUDE_DIR)

if (CppHoafParser_FOUND AND NOT TARGET CppHoafParser::CppHoafParser)
    add_library(CppHoafParser INTERFACE IMPORTED)
    set_target_properties(CppHoafParser PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${CppHoafParser_INCLUDE_DIR}/consumer;${CppHoafParser_INCLUDE_DIR}/parser;${CppHoafParser_INCLUDE_DIR}/ast;${CppHoafParser_INCLUDE_DIR}/util"
            )
endif ()
