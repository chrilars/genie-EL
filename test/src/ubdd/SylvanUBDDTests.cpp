//
// Created by mrychlicki on 7/26/21.
//

#include "ubdd/BaseUBDD.hh"
#include "SylvanTestHelper.h"
#include "UBDDTests.h"
#include <gtest/gtest.h>

using namespace fairsyn;

struct SylvanUBDDTest : public UBDDTest<SylvanUBDD>, SylvanTestHelper {
    SylvanUBDD *sylvan_base;
    SylvanUBDDTest() {
        sylvan_base = new SylvanUBDD();
    }
};

TEST_F(SylvanUBDDTest, operators) { check_operators(sylvan_base); }
TEST_F(SylvanUBDDTest, var) { check_var(sylvan_base); }
TEST_F(SylvanUBDDTest, nodeIndex) { check_nodeIndex(sylvan_base); }
TEST_F(SylvanUBDDTest, nodeSize) { check_nodeSize(sylvan_base); }
TEST_F(SylvanUBDDTest, cubes) { check_cubes(sylvan_base); }
TEST_F(SylvanUBDDTest, permute) { check_permute(sylvan_base); }
TEST_F(SylvanUBDDTest, countMinterm) { check_countMinterm(sylvan_base); }
TEST_F(SylvanUBDDTest, save_load) { check_save_load(sylvan_base); }