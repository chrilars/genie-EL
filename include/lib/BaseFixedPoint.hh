/*
* FixedPoint.hh
*
*/

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>
#include <sylvan_obj.hpp>
#include <vector>

#include "lib/Arena.hh"
#include "lib/RabinAutomaton.hh"
#include "utils/TicToc.hh"

namespace fairsyn {
    /*
    * class: FixedPoint
    *
    *
    * provides the fixed point computation for the Rabin specification
    * on finite transition systems with the edge fairness condition
    *
    */
    class BaseFixedPoint {
    public:
        /* function: printTabs
         *  used to print n number of tabs on the standard I/O during printing the results
         */
        void printTabs(const int n) {
            for (int i = 0; i < n; i++) {
                std::cout << "\t";
            }
        }

        /* function: to_dec
         *  convert a number with base "base" to a decimal number
         * position 0 is MSB
         */
        size_t to_dec(const size_t base, const std::vector<size_t> number) {
            size_t N = 0;
            for (size_t i = 1; i <= number.size(); i++) {
                N += number[i - 1] * pow(base, number.size() - i);
            }
            return N;
        }

        /* function: factorial
         *  compute the factorial
         */
        size_t factorial(const size_t n) {
            size_t ans = 1;
            for (int i = n; i > 1; i--)
                ans *= i;
            return ans;
        }

        /* function: pad_zeros
         *  pad zeros to a vector "vec1" to make its size equal to n
         */
        inline std::vector<size_t> pad_zeros(const std::vector<size_t> vec1, const size_t n) {
            std::vector<size_t> vec2 = vec1;
            for (size_t i = 0; i < n; i++) {
                vec2.push_back(0);
            }
            return vec2;
        }

        /* function: check_threshold
         *  returns true if all the elements in the vector called "vec" are below a given threshold
         */
        inline bool check_threshold(const std::vector<size_t> &vec, const size_t th) {
            for (size_t i = 0; i < vec.size(); i++) {
                if (vec[i] > th) {
                    return false;
                }
            }
            return true;
        }

        inline size_t first_index(const std::vector<size_t> &vec, const size_t th) {
            size_t i = 0;
            while (vec[i] <= th)
                i++;
            return i;
        }
    }; /* close class def */
}// namespace fairsyn
