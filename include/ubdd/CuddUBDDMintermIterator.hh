/** @file CuddUBDDMintermIterator.hh
*
*  @date 10.09.2021
*  @author Mateusz Rychlicki
*/

#pragma once

#include "UBDDMintermIterator.hh"

#include "cudd.h"
#include "cuddObj.hh"
#include "dddmp.h"

namespace genie {
    typedef BDD BDDcudd;

    /**
     * @brief Implementation of Cudd version of UBDDMintermIterator
     */
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
        CuddUBDDMintermIterator(const BDDcudd& bdd, const std::vector<size_t> &ivars);
        ~CuddUBDDMintermIterator() override;

    private:
        void begin() override;
        void next() override;
    };
}// namespace genie