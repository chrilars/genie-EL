//
// Created by mrychlicki on 7/26/21.
//

#include "BaseUBDD.hh"
#include <gtest/gtest.h>

template <class UBDD>
struct UBDDTest : testing::Test {
    void check_operators(UBDD *A) {
        EXPECT_TRUE(A->one() == A->one());
        EXPECT_TRUE(A->zero() == A->zero());
        EXPECT_FALSE(A->one() == A->zero());

        check_booleans(A->zero(), A->zero());
        check_booleans(A->zero(), A->one());
        check_booleans(A->one(), A->zero());
        check_booleans(A->one(), A->one());
    }

    void check_booleans(UBDD first, UBDD second) {
        EXPECT_EQ(~first, !first);
        EXPECT_EQ(first, !(!first));
        EXPECT_EQ(first, !(!first));
        EXPECT_TRUE((first == first));
        EXPECT_TRUE(not(first == !first));
        EXPECT_TRUE(not(first != first));
        EXPECT_TRUE((first != !first));
        EXPECT_EQ(first & second, first * second);
        EXPECT_EQ(first | second, first + second);
        EXPECT_EQ(first ^ second, (first & (!second)) | ((!first) & second));
        EXPECT_EQ(first - second, first & (!second));
        EXPECT_EQ(first.one(), first | (!first));
        EXPECT_EQ(!(first | second), (!first) & (!second));

        UBDD value_1 = first;
        UBDD value_2;
        value_2 = first;
        EXPECT_EQ(value_1, value_2);

        auto third = first;
        third &= second;
        EXPECT_EQ(first & second, third);
        third = first;
        third *= second;
        EXPECT_EQ(first * second, third);
        third = first;
        third |= second;
        EXPECT_EQ(first | second, third);
        third = first;
        third += second;
        EXPECT_EQ(first + second, third);
        third = first;
        third ^= second;
        EXPECT_EQ(first ^ second, third);
        third = first;
        third -= second;
        EXPECT_EQ(first - second, third);

        EXPECT_EQ(first < second, not(first >= second));
        EXPECT_EQ(first > second, not(first <= second));
        EXPECT_EQ(first == second, not(first != second));
        EXPECT_EQ(first <= second, first < second || first == second);
        EXPECT_EQ(first >= second, first > second || first == second);
        EXPECT_EQ(first == second, first <= second && first >= second);
        EXPECT_EQ(first != second, first < second || first > second);
    }

    void check_var(UBDD *A) {
        auto var_0 = A->var();
        auto var_1 = A->var();
        auto ivar_0 = A->var(0);
        auto ivar_1 = A->var(1);
        auto ivar_41 = A->var(41);
        auto var_42 = A->var();
        auto ivar_42 = A->var(42);

        EXPECT_EQ(var_0, ivar_0);
        EXPECT_EQ(var_1, ivar_1);
        EXPECT_EQ(var_42, ivar_42);

        EXPECT_NE(var_0, var_1);
        EXPECT_NE(var_0, var_42);
        EXPECT_NE(var_1, var_42);
        EXPECT_NE(var_42, ivar_41);
    }

    void check_nodeIndex(UBDD *A) {
        auto var_0 = A->var();
        auto var_1 = A->var();
        auto ivar_41 = A->var(41);
        auto var_42 = A->var();
        auto ivar_1 = A->var(1);

        EXPECT_EQ(var_0.nodeIndex(), 0);
        EXPECT_EQ(var_1.nodeIndex(), 1);
        EXPECT_EQ(ivar_41.nodeIndex(), 41);
        EXPECT_EQ(var_42.nodeIndex(), 42);
        EXPECT_EQ(ivar_1.nodeIndex(), 1);
    }

    void check_nodeSize(UBDD *A) {
        size_t index = 0;
        EXPECT_EQ(A->nodeSize(), index++);
        auto a_0 = A->var();
        EXPECT_EQ(A->nodeSize(), index++);

        auto var = A->var();
        EXPECT_EQ(var.nodeSize(), index);
        EXPECT_EQ(A->nodeSize(), index++);
        var.var();
        EXPECT_EQ(var.nodeSize(), index);
        EXPECT_EQ(A->nodeSize(), index++);

        auto a_420 = A->var(420);
        EXPECT_EQ(A->nodeSize(), 421);
        EXPECT_EQ(var.nodeSize(), 421);

        auto a_2137 = var.var(2137);
        EXPECT_EQ(A->nodeSize(), 2138);
        EXPECT_EQ(var.nodeSize(), 2138);
    }

    void check_cubes(UBDD *A) {
        std::vector<UBDD> variables;
        variables.push_back(A->var());
        variables.push_back(A->var());
        variables.push_back(A->var());
        variables.push_back(A->var());
        std::vector<uint8_t> values_char = {1, 1, 1, 1};
        std::vector<int> values_int = {1, 1, 1, 1};
        EXPECT_EQ(A->cube(variables, values_char), A->cube(variables, values_int));
        EXPECT_EQ(A->cube(variables), A->cube(variables, values_int));
        EXPECT_EQ(A->cube(variables), A->cube(variables, values_char));

        values_char = {1, 0, 0, 1};
        values_int = {1, 1, 1, 1};
        EXPECT_TRUE(A->cube(variables, values_char) != A->cube(variables, values_int));
        EXPECT_TRUE(A->cube(variables) == A->cube(variables, values_int));
        EXPECT_TRUE(A->cube(variables) != A->cube(variables, values_char));

        values_char = {1, 0, 0, 1};
        values_int = {1, 0, 0, 1};

        EXPECT_TRUE(A->cube(variables, values_char) == A->cube(variables, values_int));
        EXPECT_TRUE(A->cube(variables) != A->cube(variables, values_int));
        EXPECT_TRUE(A->cube(variables) != A->cube(variables, values_char));
    }

    void check_permute(UBDD *A) {
        std::vector<UBDD> variables;
        variables.push_back(A->var());
        variables.push_back(A->var());
        variables.push_back(A->var());
        variables.push_back(A->var());
        variables.push_back(A->var());
        auto cube = A->cube(variables, std::vector<int>({0, 1, 0, 0, 0}));

        EXPECT_EQ(cube, cube.permute({0, 1, 2, 3, 4}, {0, 1, 2, 3, 4}));
        EXPECT_EQ(cube, cube.permute({0, 1, 2, 3, 4}, {0, 1, 2, 4, 3}));
        EXPECT_EQ(cube, cube.permute({0, 1, 2, 3, 4}, {0, 1, 3, 2, 4}));
        EXPECT_NE(cube, cube.permute({0, 1, 2, 3, 4}, {0, 4, 2, 3, 1}));

        EXPECT_EQ(cube, cube.permute({3, 4}, {4, 3}));
        EXPECT_NE(cube, cube.permute({1, 4}, {4, 1}));
        EXPECT_EQ(cube, cube.permute({1, 2}, {2, 1}).permute({1, 2}, {2, 1}));

        EXPECT_NE(cube, cube.permute({1, 2, 3}, {2, 3, 1}));
        EXPECT_NE(cube, cube.permute({1, 2, 3}, {2, 3, 1}).permute({1, 2, 3}, {2, 3, 1}));
        EXPECT_EQ(cube, cube.permute({1, 2, 3}, {2, 3, 1}).permute({1, 2, 3}, {2, 3, 1}).permute({1, 2, 3}, {2, 3, 1}));
    }

    void check_countMinterm(UBDD *A) {
        std::vector<UBDD> variables;
        UBDD a, b, c, tautology;

        a = A->var();
        EXPECT_FLOAT_EQ(a.countMinterm(0), 0.5);
        EXPECT_FLOAT_EQ(a.countMinterm(1), 1);
        EXPECT_FLOAT_EQ(a.countMinterm(2), 2);
        EXPECT_FLOAT_EQ(a.countMinterm(3), 4);

        variables = {A->var(0), A->var(0)};
        a = A->cube(variables);
        EXPECT_FLOAT_EQ(a.countMinterm(0), 0.5);
        EXPECT_FLOAT_EQ(a.countMinterm(1), 1);
        EXPECT_FLOAT_EQ(a.countMinterm(2), 2);
        EXPECT_FLOAT_EQ(a.countMinterm(3), 4);

        variables = {A->var(0), !A->var(0)};
        a = A->cube(variables);
        EXPECT_FLOAT_EQ(a.countMinterm(0), 0);
        EXPECT_FLOAT_EQ(a.countMinterm(1), 0);
        EXPECT_FLOAT_EQ(a.countMinterm(2), 0);
        EXPECT_FLOAT_EQ(a.countMinterm(3), 0);

        variables = {A->var(0), A->var(1)};
        a = A->cube(variables, std::vector<int>({0, 1}));
        EXPECT_FLOAT_EQ(a.countMinterm(0), 0.25);
        EXPECT_FLOAT_EQ(a.countMinterm(1), 0.5);
        EXPECT_FLOAT_EQ(a.countMinterm(2), 1);
        EXPECT_FLOAT_EQ(a.countMinterm(3), 2);


        a = A->var() | A->var();
        EXPECT_FLOAT_EQ(a.countMinterm(0), 0.75);
        EXPECT_FLOAT_EQ(a.countMinterm(1), 1.5);
        EXPECT_FLOAT_EQ(a.countMinterm(2), 3);
        EXPECT_FLOAT_EQ(a.countMinterm(3), 6);

        a = A->var(42);
        b = A->var(2137);
        c = A->var();

        tautology = !((!(a & b)) ^ ((!a) | (!b)));

        EXPECT_FLOAT_EQ(tautology.countMinterm(0), 1);
        EXPECT_FLOAT_EQ((!tautology).countMinterm(0), 0);
        //        EXPECT_FLOAT_EQ(tautology.countMinterm(128), INFINITY);
        EXPECT_FLOAT_EQ((!tautology).countMinterm(128), 0);

        tautology = a | (!a) | b | c;
        EXPECT_FLOAT_EQ(tautology.countMinterm(0), 1);

        tautology = b | c | a | (!a);
        EXPECT_FLOAT_EQ(tautology.countMinterm(0), 1);


        UBDD bnf = ((!a) & b & c) | (a & (!b) & c);
        EXPECT_FLOAT_EQ(bnf.countMinterm(1), 0.5);
        EXPECT_FLOAT_EQ(bnf.countMinterm(2), 1);
        EXPECT_FLOAT_EQ(bnf.countMinterm(3), 2);
        EXPECT_FLOAT_EQ((!bnf).countMinterm(3), 6);
    }

    void check_save_load(UBDD *A) {

        UBDD a, b, c;
        check_save_and_load(A->one(), true);
        check_save_and_load(A->zero(), true);
        check_save_and_load(A->var(1));
        check_save_and_load(A->var(0));

        a = A->var(0);
        b = A->var(2);
        c = A->var(1);

        check_save_and_load(a);
        check_save_and_load(a | (!a), true);
        check_save_and_load(a & b & c);
        check_save_and_load(!((!(a & b)) ^ ((!a) | (!b))), true);
        check_save_and_load(((!(a & b)) ^ ((!a) | (!b))), true);
        check_save_and_load(((!a) & b & c) | (a & (!b) & c));
    }

    void check_save_and_load(UBDD ubdd, bool constant = false) {
        FILE *file = fopen("UBDDTest_saveAndLoad.bdd", "w");
        ubdd.save(file);
        fclose(file);

        file = fopen("UBDDTest_saveAndLoad.bdd", "r");
        UBDD fromFile = ubdd.load(file, {}, 0);
        fclose(file);
        EXPECT_EQ(ubdd, fromFile);


        file = fopen("UBDDTest_saveAndLoad.bdd", "r");
        fromFile = ubdd.load(file, {3, 4, 5}, 1);
        fclose(file);
        if (not constant) {
            EXPECT_NE(ubdd, fromFile);
        }
        EXPECT_EQ(ubdd, fromFile.permute({0, 1, 2, 3, 4, 5}, {3, 4, 5, 0, 1, 2}));
    }
};
