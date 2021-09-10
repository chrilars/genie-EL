//
// Created by mrychlicki on 7/26/21.
//

#include "UBDDTests.h"
#include "ubdd/CuddUBDD.hh"
#include <gtest/gtest.h>

using namespace fairsyn;

struct CuddUBDDTest : public UBDDTest<CuddUBDD> {
    CuddUBDD *cudd_base;
    CuddUBDDTest() {
        cudd_base = new CuddUBDD();
    }
};

TEST_F(CuddUBDDTest, CuddUBDD_operators) { check_operators(cudd_base); }
TEST_F(CuddUBDDTest, CuddUBDD_var) { check_var(cudd_base); }
TEST_F(CuddUBDDTest, CuddUBDD_nodeIndex) { check_nodeIndex(cudd_base); }
TEST_F(CuddUBDDTest, CuddUBDD_nodeSize) { check_nodeSize(cudd_base); }
TEST_F(CuddUBDDTest, CuddUBDD_cubes) { check_cubes(cudd_base); }
TEST_F(CuddUBDDTest, CuddUBDD_permute) { check_permute(cudd_base); }
TEST_F(CuddUBDDTest, CuddUBDD_countMinterm) { check_countMinterm(cudd_base); }
TEST_F(CuddUBDDTest, CuddUBDD_save_load) { check_save_load(cudd_base); }