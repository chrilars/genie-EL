message("CMAKELIST: Tests")

include(GoogleTest) # For Terminal
find_package(GTest REQUIRED) # For Clion

add_executable(FairSynTests
        ubdd/UBDDTests.h
        ubdd/UBDDMintermIteratorTests.h
        ubdd/CuddUBDDTests.cpp
        ubdd/CuddUBDDMintermIteratorTests.cpp
        ubdd/SylvanUBDDMintermIteratorTests.cpp
        ubdd/SylvanUBDDTests.cpp
        tests.cpp)

gtest_add_tests(TARGET FairSynTests)

set_target_properties(FairSynTests PROPERTIES LINKER_LANGUAGE CXX)


add_library(cudd STATIC IMPORTED)
set_property(TARGET cudd PROPERTY IMPORTED_LOCATION ${CUDD_LIBRARY_PATH}/libcudd.a)
target_link_libraries(FairSynTests cudd)

add_library(sylvan STATIC IMPORTED)
set_property(TARGET sylvan PROPERTY IMPORTED_LOCATION ${LD_LIBRARY_PATH}/libSylvan.a)
target_link_libraries(FairSynTests sylvan)

target_link_libraries(FairSynTests GTest::GTest GTest::Main)
