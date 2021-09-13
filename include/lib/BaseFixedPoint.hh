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

        /* function: check_threshold
         *  returns true if all the elements in the vector called "vec" are below a given threshold
         */
        inline bool check_threshold(const std::vector<size_t> &vec, const size_t th) {
            for (size_t i = 0; i < vec.size(); i++)
                if (vec[i] > th)
                    return false;
            return true;
        }
    }; /* close class def */
}// namespace fairsyn
