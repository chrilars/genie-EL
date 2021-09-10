//
// Created by mrychlicki on 7/26/21.
//

#include "ubdd/BaseUBDD.hh"
#include "ubdd/UBDDMintermIterator.hh"
#include <algorithm>
#include <gtest/gtest.h>
#include <vector>

using namespace fairsyn;

template <class UBDD>
struct UBDDMintermIteratorTests : testing::Test {
    UBDD base;
    void check_done() {
        simple_test_done(base.var(0), {0}, 1);
        simple_test_done(!base.var(0), {0}, 1);
        simple_test_done(base.var(0) | base.var(1), {0, 1}, 3);
        simple_test_done(base.var(0) | (!base.var(0)), {0}, 2);
        simple_test_done(base.var(0) & (!base.var(0)), {0}, 0);

        simple_test_done(base.var(42), {42}, 1);
        simple_test_done(!base.var(42), {42}, 1);
        simple_test_done(base.var(42) | base.var(2137), {42, 2137}, 3);
        simple_test_done(base.var(42) | (!base.var(42)), {42}, 2);
        simple_test_done(base.var(42) & (!base.var(42)), {42}, 0);

        UBDD random = (base.var(42) ^ (!base.var(2137))) | base.var(69);
        simple_test_done(random, {69, 42, 2137}, 6);
        simple_test_done(!random, {69, 42, 2137}, 2);
    }

    void check_numberOfMinterms() {
        simple_test_numberOfMinterms(base.var(0), {0}, 1);
        simple_test_numberOfMinterms(!base.var(0), {0}, 1);
        simple_test_numberOfMinterms(base.var(0) | base.var(1), {0, 1}, 3);
        simple_test_numberOfMinterms(base.var(0) | (!base.var(0)), {0}, 2);
        simple_test_numberOfMinterms(base.var(0) & (!base.var(0)), {0}, 0);

        simple_test_numberOfMinterms(base.var(42), {42}, 1);
        simple_test_numberOfMinterms(!base.var(42), {42}, 1);
        simple_test_numberOfMinterms(base.var(42) | base.var(2137), {42, 2137}, 3);
        simple_test_numberOfMinterms(base.var(42) | (!base.var(42)), {42}, 2);
        simple_test_numberOfMinterms(base.var(42) & (!base.var(42)), {42}, 0);

        UBDD random = (base.var(42) ^ (!base.var(2137))) | base.var(69);
        simple_test_numberOfMinterms(random, {69, 42, 2137}, 6);
        simple_test_done(!random, {69, 42, 2137}, 2);

        UBDD large_ubdd = base.var(0);
        std::vector<size_t> ivars_large(100);
        for (int i = 0; i < 100; i++) {
            large_ubdd |= base.var(20 * i);
            ivars_large[i] = 20 * i;
        }
        EXPECT_ANY_THROW(large_ubdd.generateMintermIterator(ivars_large));
    }

    void check_currentMinterm() {
        simple_test_currentMinterm(base.var(0), {0}, {{1}});
        simple_test_currentMinterm(!base.var(0), {0}, {{0}});
        simple_test_currentMinterm(base.var(0) | base.var(1), {0, 1}, {
                                                                              {1, 1},
                                                                              {0, 1},
                                                                              {1, 0},
                                                                      });
        simple_test_currentMinterm(base.var(0) | (!base.var(0)), {0}, {{1}, {0}});
        simple_test_currentMinterm(base.var(0) & (!base.var(0)), {0}, {});

        simple_test_currentMinterm(base.var(42), {42}, {{1}});
        simple_test_currentMinterm(!base.var(42), {42}, {{0}});
        simple_test_currentMinterm(base.var(42) | base.var(2137), {42, 2137}, {
                                                                                      {1, 1},
                                                                                      {0, 1},
                                                                                      {1, 0},
                                                                              });
        simple_test_currentMinterm(base.var(42) | (!base.var(42)), {42}, {{1}, {0}});
        simple_test_currentMinterm(base.var(42) & (!base.var(42)), {42}, {});

        UBDD random = (base.var(42) ^ (!base.var(2137))) | base.var(69);
        simple_test_currentMinterm(random, {69, 42, 2137}, {
                                                                   {1, 1, 1},
                                                                   {1, 1, 0},
                                                                   {1, 0, 1},
                                                                   {1, 0, 0},
                                                                   {0, 1, 1},
                                                                   {0, 0, 0},
                                                           });
        simple_test_currentMinterm(!random, {69, 42, 2137}, {
                                                                    {0, 1, 0},
                                                                    {0, 0, 1},
                                                            });
    }

private:
    void simple_test_done(UBDD cube, std::vector<size_t> ivars, size_t prediction) {
        UBDDMintermIterator *it = cube.generateMintermIterator(ivars);
        for (size_t i = 0; i < prediction; i++) {
            EXPECT_FALSE(it->done());
            it->operator++();
        }
        EXPECT_TRUE(it->done());
    }

    void simple_test_numberOfMinterms(UBDD cube, std::vector<size_t> ivars, double prediction) {
        UBDDMintermIterator *it = cube.generateMintermIterator(ivars);
        EXPECT_FLOAT_EQ(it->numberOfMinterms(), prediction);
    }

    void simple_test_currentMinterm(UBDD cube, std::vector<size_t> ivars, std::set<std::vector<size_t>> prediction) {
        UBDDMintermIterator *it = cube.generateMintermIterator(ivars);
        std::set<std::vector<size_t>> s;
        std::vector<size_t> minterm;
        for (size_t i = 0; i < prediction.size(); i++) {
            EXPECT_FALSE(it->done());
            minterm = it->shortMinterm();
            auto pointer = prediction.find(minterm);
            EXPECT_FALSE(pointer == prediction.end());
            EXPECT_TRUE(s.find(minterm) == s.end());
            s.insert(minterm);
            it->operator++();
        }
        EXPECT_TRUE(it->done());
    }
};
