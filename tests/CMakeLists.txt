message("CMAKELIST: FaiSyn.tests")

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

target_link_libraries(FairSynTests GTest::GTest GTest::Main)
target_link_libraries(FairSynTests UBDD)
