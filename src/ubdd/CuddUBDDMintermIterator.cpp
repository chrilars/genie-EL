/** @file CuddUBDDMintermIterator.cpp
*
*  @date 10.09.2021
*  @author Mateusz Rychlicki
*/

#include "ubdd/CuddUBDDMintermIterator.hh"

namespace genie {
    CuddUBDDMintermIterator::CuddUBDDMintermIterator(const BDDcudd &bdd, const std::vector<size_t> &ivars) {
        ivars_.assign(ivars.begin(), ivars.end());
        nvars_ = ivars.size();
        done_ = false;
        counter_ = 0;
        max_ivars_ = *max_element(ivars_.begin(), ivars_.end());
        minterm_ = std::vector<size_t>(max_ivars_ + 1, 2);

        if (nvars_ > sizeof(size_t) * 8) {
            throw std::invalid_argument(
                    "Error: UBDDMintermIterator: number of variables we iterate over is limited to highest number in size_t.");
        }
        bdd_ = bdd;
        /* check if bdd depends only on bdd variables with indices ivars */
        std::vector<unsigned int> sup = bdd.SupportIndices();
        for (size_t i = 0; i < sup.size(); i++) {
            int marker = 0;
            for (size_t j = 0; j < nvars_; j++) {
                if (sup[i] == ivars_[j]) {
                    marker = 1;
                    break;
                }
            }
            if (!marker) {
                throw std::invalid_argument(
                        "Error: CuddUBDDMintermIterator: the bdd depends on variables with index outside ivars.");
            }
        }

        begin();
    }

    CuddUBDDMintermIterator::~CuddUBDDMintermIterator() {
        Cudd_GenFree(gen_);
    }

    void CuddUBDDMintermIterator::begin() {
        /* initialize gen_ and grab first cube */
        DdManager *manager = bdd_.manager();
        DdNode *node = bdd_.getNode();
        gen_ = Cudd_FirstCube(manager, node, &cube_, &value_);
        if (Cudd_IsGenEmpty(gen_)) {
            done_ = true;
            return;
        }
        /* determine the number of minterms */
        nminterm_ = bdd_.CountMinterm(nvars_);
        progress_ = 1;

        /* check cube for don't cares */
        idontcares_.reserve(nvars_);
        ndontcares_ = 0;
        nexpand_ = 1;
        iterator_ = 0;
        /* indices of the variables with don't cares */
        for (size_t i = 0; i < nvars_; i++) {
            if (cube_[ivars_[i]] == 2) {
                idontcares_[ndontcares_] = ivars_[i];
                ndontcares_++;
                nexpand_ *= 2;
            }
        }
        /* if the cube contins don't cares, start expanding the cube to a minterm */
        if (ndontcares_) {
            /* start with all zero entries */
            for (size_t i = 0; i < ndontcares_; i++)
                cube_[idontcares_[i]] = 0;
        }

        for (unsigned long ivar: ivars_)
            minterm_[ivar] = cube_[ivar];
    }

    void CuddUBDDMintermIterator::next() {
        iterator_++;
        progress_++;
        /* get new cube or expand cube further */
        if (iterator_ < nexpand_) {
            /* still not don with this cube */
            /* generate binary value for iterator_ to fill cube */
            std::bitset<sizeof(size_t) * 8> iterator2binary(iterator_);
            for (size_t i = 0; i < ndontcares_; i++)
                cube_[idontcares_[i]] = iterator2binary[i];
        } else {
            /* restore the cube (put the don't cares) */
            for (size_t i = 0; i < ndontcares_; i++)
                cube_[idontcares_[i]] = 2;
            /* now get the next cube */
            Cudd_NextCube(gen_, &cube_, &value_);
            if (Cudd_IsGenEmpty(gen_)) {
                done_ = true;
                return;
            }
            /* check for don't cares */
            ndontcares_ = 0;
            nexpand_ = 1;
            iterator_ = 0;
            for (size_t i = 0; i < nvars_; i++) {
                if (cube_[ivars_[i]] == 2) {
                    idontcares_[ndontcares_] = ivars_[i];
                    ndontcares_++;
                    nexpand_ *= 2;
                }
            }
            /* if the cube contins don't cares, start expanding the cube to a minterm */
            if (ndontcares_) {
                /* start with all zero entries */
                for (size_t i = 0; i < ndontcares_; i++)
                    cube_[idontcares_[i]] = 0;
            }
        }
        for (unsigned long ivar: ivars_)
            minterm_[ivar] = cube_[ivar];
    }
}// namespace genie
