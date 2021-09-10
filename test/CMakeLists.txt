include(GoogleTest) # For Terminal
find_package(GTest REQUIRED) # For Clion

# Now simply link against gtest as needed. Eg
message(mascot/fairsyn/test)
add_executable(fairsyn_test
        tests.cpp
        src/ubdd/CuddUBDDTests.cpp
        src/ubdd/CuddUBDDMintermIteratorTests.cpp
        src/ubdd/SylvanTestHelper.h
        src/ubdd/SylvanUBDDTests.cpp
        src/ubdd/SylvanUBDDMintermIteratorTests.cpp
        src/ubdd/UBDDTests.h
        src/ubdd/UBDDMintermIteratorTests.h
        )
target_compile_features(fairsyn_test PRIVATE cxx_auto_type)
target_link_libraries(fairsyn_test FairSyn::fairsyn GTest::GTest GTest::Main)
