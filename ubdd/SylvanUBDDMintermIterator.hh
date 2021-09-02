/*
 * SylvanUBDDMintermIterator.hh
 *
 *  created on: 12.07.2021
 *      author: Mateusz Rychlicki
 *
 */

#ifndef SYLVANUBDDMINTERMITERATOR_HH_
#define SYLVANUBDDMINTERMITERATOR_HH_

#include "UBDDMintermIterator.hh"
#include "sylvan_obj.hpp"
#include <vector>

#include "sylvan_bdd.h"

typedef sylvan::Bdd BDDsylvan;

class SylvanUBDDMintermIterator : public UBDDMintermIterator {
    sylvan::MTBDD leaf_;
    sylvan::MTBDD dd_;
    sylvan::MTBDD variables_;
    std::vector<size_t> sorted_ivars_;
    uint8_t *arr_;

public:
    SylvanUBDDMintermIterator(sylvan::MTBDD dd, sylvan::MTBDD variables, const std::vector<size_t> &ivars) {
        dd_ = dd;
        variables_ = variables;
        sylvan::sylvan_protect(&dd_);
        sylvan::sylvan_protect(&variables_);
        ivars_ = ivars;
        sorted_ivars_ = ivars;
        sort(sorted_ivars_.begin(), sorted_ivars_.end());
        nvars_ = ivars.size();
        done_ = false;
        counter_ = 0;
        arr_ = new uint8_t[nvars_];
        max_ivars_ = *max_element(ivars_.begin(), ivars_.end() - 1);
        minterm_ = std::vector<size_t>(max_ivars_ + 1, 2);

        if (nvars_ > sizeof(size_t) * 8) {
            std::ostringstream os;
            os << "Error: UBDDMintermIterator: number of variables we iterate over is limited to highest number in size_t.";
            throw std::invalid_argument(os.str().c_str());
        }

        begin();
    }

    ~SylvanUBDDMintermIterator() {
        sylvan::sylvan_unprotect(&dd_);
        sylvan::sylvan_unprotect(&variables_);
    }

protected:
    void begin() {
        nminterm_ = BDDsylvan(dd_).SatCount(sylvan::BddSet(sylvan::Bdd(variables_)));
        leaf_ = sylvan::mtbdd_enum_all_first(dd_, variables_, arr_, NULL);
        done_ = (leaf_ == sylvan::mtbdd_false);
        for (size_t i = 0; i < ivars_.size(); i++)
            minterm_[sorted_ivars_[i]] = arr_[i];
    }

    void next() {
        leaf_ = sylvan::mtbdd_enum_all_next(dd_, variables_, arr_, NULL);
        done_ = (leaf_ == sylvan::mtbdd_false);
        for (size_t i = 0; i < ivars_.size(); i++)
            minterm_[sorted_ivars_[i]] = arr_[i];
    }
};


#endif /* SYLVANUBDDMINTERMITERATOR_HH_ */
