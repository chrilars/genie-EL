/*
 * SylvanUBDD.hh
 *
 *  created on: 01.07.2021
 *      author: Mateusz Rychlicki
 *
 */
#pragma once

#include "BaseUBDD.hh"
#include "SylvanUBDDMintermIterator.hh"
#include <cfloat>
#include <cmath>
#include <map>
#include <memory>


using namespace sylvan;
namespace fairsyn {
    class SylvanUBDD : public BaseUBDD<SylvanUBDD> {

    public:
        inline static size_t size_ = 0;
        inline static std::map<size_t, SylvanUBDD *> nodes_map = std::map<size_t, SylvanUBDD *>();
        BDDsylvan bdd_;
        SylvanUBDD();
        SylvanUBDD(BDDsylvan const &from);
        SylvanUBDD(SylvanUBDD const &from);
        SylvanUBDD zero() const override;
        SylvanUBDD one() const override;
        SylvanUBDD var() const override;
        SylvanUBDD var(size_t index) const override;
        size_t nodeSize() const override;
        size_t nodeIndex() const override;
        SylvanUBDD cube(std::vector<SylvanUBDD> variables, std::vector<uint8_t> values) const override;
        SylvanUBDD cube(std::vector<SylvanUBDD> variables, std::vector<int> values) const override;
        SylvanUBDD cube(std::vector<SylvanUBDD> variables) const override;
        SylvanUBDD permute(const std::vector<size_t> &from, const std::vector<size_t> &to) const override;
        double countMinterm(size_t nvars) const override;
        void printMinterm() const override;
        SylvanUBDD existAbstract(const SylvanUBDD &cube) const override;
        SylvanUBDD andAbstract(const SylvanUBDD &g, const SylvanUBDD &cube) const override;
        UBDDMintermIterator *generateMintermIterator(std::vector<size_t> &ivars) const override;
        SylvanUBDD complement(SylvanUBDD &symbolicSet,
                              std::vector<size_t> nofGridPoints,
                              std::vector<size_t> nofBddVars,
                              std::vector<std::vector<size_t>> indBddVars,
                              size_t dim) override;
        SylvanUBDD computePolytope(const size_t p,
                                   const std::vector<double> &H,
                                   const std::vector<double> &h,
                                   int type,
                                   size_t dim,
                                   const std::vector<double> &eta,
                                   const std::vector<double> &z,
                                   const std::vector<double> &firstGridPoint,
                                   const std::vector<size_t> &nofBddVars,
                                   const std::vector<std::vector<size_t>> &indBddVars,
                                   const std::vector<size_t> &nofGridPoints) override;
        int save(FILE *file) override;
        SylvanUBDD load(FILE *file, std::vector<int> composeids, int newID) override;
        SylvanUBDD transfer(const SylvanUBDD &destination) const override;
        bool isCoverEqual(const SylvanUBDD &other) const override;
        SylvanUBDD &operator=(const SylvanUBDD &right);
        bool operator==(const SylvanUBDD &other) const override;
        bool operator!=(const SylvanUBDD &other) const override;
        bool operator<=(const SylvanUBDD &other) const override;
        bool operator>=(const SylvanUBDD &other) const override;
        bool operator<(const SylvanUBDD &other) const override;
        bool operator>(const SylvanUBDD &other) const override;
        SylvanUBDD operator!() const override;
        SylvanUBDD operator~() const override;
        SylvanUBDD operator*(const SylvanUBDD &other) const override;
        SylvanUBDD operator*=(const SylvanUBDD &other) override;
        SylvanUBDD operator&(const SylvanUBDD &other) const override;
        SylvanUBDD operator&=(const SylvanUBDD &other) override;
        SylvanUBDD operator+(const SylvanUBDD &other) const override;
        SylvanUBDD operator+=(const SylvanUBDD &other) override;
        SylvanUBDD operator|(const SylvanUBDD &other) const override;
        SylvanUBDD operator|=(const SylvanUBDD &other) override;
        SylvanUBDD operator^(const SylvanUBDD &other) const override;
        SylvanUBDD operator^=(const SylvanUBDD &other) override;
        SylvanUBDD operator-(const SylvanUBDD &other) const override;
        SylvanUBDD operator-=(const SylvanUBDD &other) override;

    private:
        static std::vector<BDDsylvan> extractBDDs(const std::vector<SylvanUBDD> &ubdds) {
            std::vector<BDDsylvan> bdds(ubdds.size());
            for (size_t i = 0; i < ubdds.size(); i++)
                bdds[i] = ubdds[i].bdd_;
            return bdds;
        }
    };
}// namespace fairsyn