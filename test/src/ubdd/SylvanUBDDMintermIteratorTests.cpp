//
// Created by mrychlicki on 7/26/21.
//

#include "SylvanTestHelper.h"
//#include "SylvanUBDD.hh"
#include "UBDDMintermIteratorTests.h"

using namespace fairsyn;

struct SylvanUBDDMintermIteratorTests : UBDDMintermIteratorTests<SylvanUBDD>, SylvanTestHelper {
};

TEST_F(SylvanUBDDMintermIteratorTests, done) { check_done(); }
TEST_F(SylvanUBDDMintermIteratorTests, numberOfMinterms) { check_numberOfMinterms(); }
TEST_F(SylvanUBDDMintermIteratorTests, currentMinterm) { check_currentMinterm(); }
