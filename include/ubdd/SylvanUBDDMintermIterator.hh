/*
 * SylvanUBDDMintermIterator.hh
 *
 *  created on: 12.07.2021
 *      author: Mateusz Rychlicki
 *
 */

#pragma once

#include "UBDDMintermIterator.hh"
#include "sylvan_obj.hpp"
#include <vector>

#include "sylvan_bdd.h"

namespace fairsyn {
    typedef sylvan::Bdd BDDsylvan;

    class SylvanUBDDMintermIterator : public UBDDMintermIterator {
        sylvan::MTBDD leaf_;
        sylvan::MTBDD dd_;
        sylvan::MTBDD variables_;
        std::vector<size_t> sorted_ivars_;
        uint8_t *arr_;

    public:
        SylvanUBDDMintermIterator(sylvan::MTBDD dd, sylvan::MTBDD variables, const std::vector<size_t> &ivars);
        ~SylvanUBDDMintermIterator();

    protected:
        void begin();
        void next();
    };
}// namespace fairsyn
