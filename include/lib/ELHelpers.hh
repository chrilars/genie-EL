#pragma once

#include <algorithm>
#include <cassert>
#include <iostream>
#include <cmath>
#include <iterator>
#include <vector>
#include "structures.hh"

namespace ELHelpers {

    inline bool proper_subset(const std::vector<bool>& l1, const std::vector<bool>& l2) {
        // Check if l1 subset of l2
        // forall i in l1: if l1 then also l2
        // but there must exist at least one index in l2 which is not in l1
        bool flag = false;
        for (size_t i = 0; i < l1.size(); ++i) {
            if (l2[i] && !l1[i])
                flag = true;
            if (l1[i] && !l2[i])
                return false;
        }
        return flag;
    }

    inline void print_label(ZielonkaNode *z) {
        for (size_t i = 0; i < z->label.size(); ++i){
            if (z->label[i]) {
                std::cout << static_cast<char>('a' + i) << ", ";
            } else {
                std::cout << "   ";
            }
        } //std::cout << '\n';
    }


    inline std::vector<bool> label_difference(const std::vector<bool>& t, const std::vector<bool>& s) {
        // compute difference t-s, where s,t are vectors of colors
        size_t size = t.size();
        std::vector<bool> difference(size, false);
        for (size_t i = 0; i < size; ++i)
            difference[i] = t[i] && !s[i];
        return difference;
    }


    // Standard algorithm for powerset generation: https://www.geeksforgeeks.org/power-set/
    inline std::vector<std::vector<bool>> powerset(size_t colors) {
        // remake this so that it works with new evaluation function (std::vector<bool>)
        size_t ps_size = pow(2, colors);
        std::vector<std::vector<bool>> ps;
        for (size_t count = 0; count < ps_size; ++count) {
            std::vector<bool> tmp(colors, 0);
            for (size_t bit = 0; bit < colors; ++bit) {
                if (count & (1 << bit))
                    tmp[bit] = true;
            }
            ps.push_back(tmp);
        }
        return ps;
    }

    inline std::vector<size_t> preprocess_to_UBDD(const std::vector<bool>& label) {
        std::vector<size_t> preprocessed;
        for (size_t i = 0; i < label.size(); ++i){
            if (label[i]) preprocessed.push_back(i);
        } return preprocessed;
    }
} // namespace EMHelpers
