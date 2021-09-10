/*
* CuddUBDDMintermIterator.hh
*
*  created on: 11.07.2021
*      author: Mateusz Rychlicki
*/
#pragma once

#include "UBDDMintermIterator.hh"

#include "cudd.h"
#include "cuddObj.hh"
#include "dddmp.h"
namespace fairsyn {
    typedef BDD BDDcudd;

    class CuddUBDDMintermIterator : public UBDDMintermIterator {
    private:
        BDDcudd bdd_;
        int *cube_;
        DdGen *gen_;
        CUDD_VALUE_TYPE value_;
        std::vector<size_t> idontcares_;
        size_t ndontcares_;
        size_t iterator_;
        size_t nexpand_;

    public:
        CuddUBDDMintermIterator(BDDcudd bdd, const std::vector<size_t> &ivars);
        ~CuddUBDDMintermIterator();

    private:
        void begin();
        void next();
    };
}// namespace fairsyn