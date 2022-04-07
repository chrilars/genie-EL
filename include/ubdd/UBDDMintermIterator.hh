/** @file UBDDMintermIterator.hh
*
*  @date 10.09.2021
*  @author Mateusz Rychlicki
*/

#pragma once

#include <algorithm>
#include <bitset>
#include <cassert>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace fairsyn {
    /**
     * @brief An interface for Universal BDD Minterm Iterator
     * @details
     * Miniterm iterator help with iteration over minterms
     */
    class UBDDMintermIterator {
    protected:
        std::vector<size_t> ivars_;
        std::vector<size_t> minterm_;
        bool done_;
        size_t max_ivars_;
        size_t nvars_;
        size_t progress_;
        double nminterm_;
        size_t counter_;

    public:
        virtual ~UBDDMintermIterator(){};
        /**
         * @brief Calculate next minterm.
         */
        void operator++() {
            if (!done_) {
                next();
            }
        }

        /**
         * @brief check if there is no next minterm.
         * @return %True, if there is no next minterm, otherwise %False.
         */
        inline bool done() { return done_; }

        /**
         * @brief Returns current minterm.
         * @details
         * currentMinterm.size() == max(ivars) + 1.
         * currentMinterm[i] == 0 if i \in ivars and variable i is false
         * currentMinterm[i] == 1 if i \in ivars and variable i is true
         * currentMinterm[i] == 2 if i \notin ivars 'dont care'
         *
         * @return current minterm.
         */
        inline const std::vector<size_t> &currentMinterm() {
            return minterm_;
        }

        /**
         * @brief Same as %currentMinterm but representation is based on ivars.
         * @details result[i] := currentMinterm[ivars[i]]
         *
         * @return current minterm
         */
        inline std::vector<size_t> shortMinterm() {
            std::vector<size_t> result(nvars_);
            for (size_t i = 0; i < nvars_; i++)
                result[i] = minterm_[ivars_[i]];
            return result;
        }

        /**
         * @brief Prints shortMinterm.
         */
        inline void printMinterm() {
            for (size_t i = 0; i < nvars_; i++)
                std::cout << minterm_[ivars_[i]];
            std::cout << std::endl;
        }

        /**
         * @brief return number of minterms.
         * @return number of minterms.
         */
        inline double numberOfMinterms() const {
            return nminterm_;
        }

        /**
         * @brief Prints progress.
         */
        inline void printProgress(void) {
            if ((size_t)(progress_ / nminterm_ * 100) >= counter_) {
                if ((counter_ % 10) == 0)
                    std::cout << counter_;
                else
                    std::cout << ".";
                counter_++;
            }
            std::flush(std::cout);
        }

    protected:
        /**
         * @brief initialize and calculate first minterm
         */
        virtual void begin() = 0;

        /**
         * @brief calculate next minterm
         */
        virtual void next() = 0;
    };
}// namespace fairsyn