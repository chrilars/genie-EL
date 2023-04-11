/** @file SylvanUBDDMintermIterator.hh
*
*  @date 10.09.2021
*  @author Mateusz Rychlicki
*/

#pragma once

#include "UBDDMintermIterator.hh"
#include "sylvan_obj.hpp"
#include <vector>

#include "sylvan_bdd.h"

namespace genie {
    typedef sylvan::Bdd BDDsylvan;

    /**
     * @brief Implementation of Sylvan version of UBDDMintermIterator
     */
    class SylvanUBDDMintermIterator : public UBDDMintermIterator {
        sylvan::MTBDD leaf_;
        sylvan::MTBDD dd_;
        sylvan::MTBDD variables_;
        std::vector<size_t> sorted_ivars_;
        uint8_t *arr_;

    public:
        SylvanUBDDMintermIterator(sylvan::MTBDD dd, sylvan::MTBDD variables, const std::vector<size_t> &ivars);
        ~SylvanUBDDMintermIterator() override;

    protected:
        void begin() override;
        void next() override;
    };
}// namespace genie
