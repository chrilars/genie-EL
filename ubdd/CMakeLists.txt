cmake_minimum_required(VERSION 3.9.1)
message("CMAKELIST: FaiSyn.UBDD")
project(UBDD VERSION 1.0.1 DESCRIPTION "Universal Binary Decision Diagrams ")
add_library(UBDD STATIC # STATIC OR SHARED
        BaseUBDD.hh
        CuddUBDD.hh
        SylvanUBDD.hh
        UBDDMintermIterator.hh
        CuddUBDDMintermIterator.hh
        SylvanUBDDMintermIterator.hh)

set_target_properties(UBDD PROPERTIES LINKER_LANGUAGE CXX VERSION ${PROJECT_VERSION})

add_library(cudd STATIC IMPORTED)
set_property(TARGET cudd PROPERTY IMPORTED_LOCATION ${CUDD_LIBRARY_PATH}/libcudd.a)
target_link_libraries(UBDD cudd)

add_library(sylvan STATIC IMPORTED)
set_property(TARGET sylvan PROPERTY IMPORTED_LOCATION ${SYLVAN_PATH}/libSylvan.a)
target_link_libraries(UBDD sylvan)

#target_include_directories(UBDD PRIVATE .)