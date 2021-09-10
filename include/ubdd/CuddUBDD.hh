/*
 * CuddUBDD.hh
 *
 *  created on: 01.07.2021 (In Progress)
 *      author: Mateusz Rychlicki
 *
 */
#pragma once

#include "BaseUBDD.hh"
#include "CuddUBDDMintermIterator.hh"
#include <complex>
#include <memory>
#include <utility>

#include "cudd.h"
#include "cuddObj.hh"
#include "dddmp.h"

namespace fairsyn {
    class CuddUBDD : public BaseUBDD<CuddUBDD> {
    public:
        BDDcudd bdd_;
        std::shared_ptr<Cudd> cudd_;

        CuddUBDD(BDDcudd &bdd, std::shared_ptr<Cudd> cudd);
        CuddUBDD(unsigned int numVars,
                 unsigned int numVarsZ,
                 unsigned int numSlots = CUDD_UNIQUE_SLOTS,
                 unsigned int cacheSize = CUDD_CACHE_SLOTS,
                 unsigned long maxMemory = 0,
                 PFC defaultHandler = defaultError);
        CuddUBDD(const Cudd &x);
        CuddUBDD(CuddUBDD &ubdd);
        CuddUBDD(const CuddUBDD &ubdd);
        CuddUBDD();
        CuddUBDD zero() const override;
        CuddUBDD one() const override;
        CuddUBDD var() const override;
        CuddUBDD var(size_t index) const override;
        size_t nodeSize() const override;
        size_t nodeIndex() const override;
        size_t nodeCount() const;
        CuddUBDD cube(std::vector<CuddUBDD> variables, std::vector<int> values) const override;
        CuddUBDD cube(std::vector<CuddUBDD> variables, std::vector<uint8_t> values) const override;
        CuddUBDD cube(std::vector<CuddUBDD> variables) const override;
        CuddUBDD permute(const std::vector<size_t> &from, const std::vector<size_t> &to) const override;
        double countMinterm(size_t nvars) const override;
        void printMinterm() const override;
        double readEpsilon() const;
        CuddUBDD existAbstract(const CuddUBDD &cube) const override;
        CuddUBDD andAbstract(const CuddUBDD &g, const CuddUBDD &cube) const override;
        UBDDMintermIterator *generateMintermIterator(std::vector<size_t> &ivars) const override;
        CuddUBDD complement(CuddUBDD &symbolicSet,
                            std::vector<size_t> nofGridPoints,
                            std::vector<size_t> nofBddVars,
                            std::vector<std::vector<size_t>> indBddVars,
                            size_t dim) override;
        CuddUBDD computePolytope(const size_t p,
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
        CuddUBDD load(FILE *file, std::vector<int> composeids, int newID) override;
        CuddUBDD transfer(const CuddUBDD &destination) const override;
        bool isCoverEqual(const CuddUBDD &other) const override;
        CuddUBDD &operator=(const CuddUBDD &right) override;
        bool operator==(const CuddUBDD &other) const override;
        bool operator!=(const CuddUBDD &other) const override;
        bool operator<=(const CuddUBDD &other) const override;
        bool operator>=(const CuddUBDD &other) const override;
        bool operator<(const CuddUBDD &other) const override;
        bool operator>(const CuddUBDD &other) const override;
        CuddUBDD operator!() const override;
        CuddUBDD operator~() const override;
        CuddUBDD operator*(const CuddUBDD &other) const override;
        CuddUBDD operator*=(const CuddUBDD &other) override;
        CuddUBDD operator&(const CuddUBDD &other) const override;
        CuddUBDD operator&=(const CuddUBDD &other) override;
        CuddUBDD operator+(const CuddUBDD &other) const override;
        CuddUBDD operator+=(const CuddUBDD &other) override;
        CuddUBDD operator|(const CuddUBDD &other) const override;
        CuddUBDD operator|=(const CuddUBDD &other) override;
        CuddUBDD operator^(const CuddUBDD &other) const override;
        CuddUBDD operator^=(const CuddUBDD &other) override;
        CuddUBDD operator-(const CuddUBDD &other) const override;
        CuddUBDD operator-=(const CuddUBDD &other) override;
        double computeEllipsoidRadius(const std::vector<double> L, size_t dim, const std::vector<double> eta, const std::vector<double> z) const;

        static std::vector<BDDcudd> extractBDDs(std::vector<CuddUBDD> ubdds) {
            std::vector<BDDcudd> bdds(ubdds.size());
            for (size_t i = 0; i < ubdds.size(); i++)
                bdds[i] = ubdds[i].bdd_;
            return bdds;
        }
    };
}// namespace fairsyn