/** @file SylvanUBDDMintermIterator.cpp
*
*  @date 10.09.2021
*  @author Mateusz Rychlicki
*/

#include "ubdd/SylvanUBDDMintermIterator.hh"

namespace genie {
    SylvanUBDDMintermIterator::SylvanUBDDMintermIterator(sylvan::MTBDD dd, sylvan::MTBDD variables, const std::vector<size_t> &ivars) {
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
        max_ivars_ = *max_element(ivars_.begin(), ivars_.end());
        minterm_ = std::vector<size_t>(max_ivars_ + 1, 2);

        if (nvars_ > sizeof(size_t) * 8) {
            throw std::invalid_argument(
                    "Error: UBDDMintermIterator: number of variables we iterate over is limited to highest number in size_t.");
        }

        begin();
    }

    SylvanUBDDMintermIterator::~SylvanUBDDMintermIterator() {
        sylvan::sylvan_unprotect(&dd_);
        sylvan::sylvan_unprotect(&variables_);
    }

    void SylvanUBDDMintermIterator::begin() {
        nminterm_ = BDDsylvan(dd_).SatCount(sylvan::BddSet(sylvan::Bdd(variables_)));
        leaf_ = sylvan::mtbdd_enum_all_first(dd_, variables_, arr_, NULL);
        done_ = (leaf_ == sylvan::mtbdd_false);
        for (size_t i = 0; i < ivars_.size(); i++)
            minterm_[sorted_ivars_[i]] = arr_[i];
        /* added by kaushik */
        progress_ = 1;
    }

    void SylvanUBDDMintermIterator::next() {
        /* added by kaushik */
        progress_++;
        /***/
        leaf_ = sylvan::mtbdd_enum_all_next(dd_, variables_, arr_, NULL);
        done_ = (leaf_ == sylvan::mtbdd_false);
        for (size_t i = 0; i < ivars_.size(); i++)
            minterm_[sorted_ivars_[i]] = arr_[i];
    }
}// namespace genie
