//
// Created by mrychlicki on 7/26/21.
//

#include "ubdd/CuddUBDD.hh"
#include "UBDDMintermIteratorTests.h"
#include <gtest/gtest.h>

using namespace fairsyn;

struct CuddUBDDMintermIteratorTests : UBDDMintermIteratorTests<CuddUBDD> {
};

TEST_F(CuddUBDDMintermIteratorTests, done) { check_done(); }
TEST_F(CuddUBDDMintermIteratorTests, numberOfMinterms) { check_numberOfMinterms(); }
TEST_F(CuddUBDDMintermIteratorTests, currentMinterm) { check_currentMinterm(); }