/*
 * SylvanUBDD.hh
 *
 *  created on: 01.07.2021
 *      author: Mateusz Rychlicki
 *
 */
#include "BaseUBDD.hh"
#include "SylvanUBDDMintermIterator.hh"
#include <cfloat>
#include <cmath>
#include <map>
#include <memory>

#ifndef SYLVANUBDD_HH_
#    define SYLVANUBDD_HH_
typedef sylvan::Bdd BDDsylvan;

using namespace sylvan;

class SylvanUBDD : public BaseUBDD<SylvanUBDD> {

public:
    inline static size_t size_ = 0;
    inline static std::map<size_t, SylvanUBDD *> nodes_map = std::map<size_t, SylvanUBDD *>();
    BDDsylvan bdd_;

    SylvanUBDD() {
    }

    SylvanUBDD(BDDsylvan const &from) {
        bdd_ = BDDsylvan(from);
    }

    SylvanUBDD(SylvanUBDD const &from) {
        bdd_ = BDDsylvan(from.bdd_);
    }

    SylvanUBDD zero() const override {
        return SylvanUBDD(bdd_.bddZero());
    }

    SylvanUBDD one() const override {
        return SylvanUBDD(bdd_.bddOne());
    }

    SylvanUBDD var() const override {
        size_t new_index = size_++;
        nodes_map[new_index] = new SylvanUBDD(bdd_.bddVar(new_index));
        return *nodes_map[new_index];
    }

    SylvanUBDD var(size_t index) const override {
        size_t new_index;
        while (size_ <= index) {
            new_index = size_++;
            nodes_map[new_index] = new SylvanUBDD(bdd_.bddVar(new_index));
        }
        return *nodes_map[index];
    }

    size_t nodeSize() const override {
        return size_;
    }

    size_t nodeIndex() const override {
        return mtbdd_getvar(bdd_.GetBDD());
    }

    SylvanUBDD cube(std::vector<SylvanUBDD> variables, std::vector<uint8_t> values) const override {
        sylvan::BddSet bddSet = BddSet::fromVector(extractBDDs(variables));
        BDDsylvan bdd = bdd_.bddCube(bddSet, values);
        return SylvanUBDD(bdd);
    }

    SylvanUBDD cube(std::vector<SylvanUBDD> variables, std::vector<int> values) const override {
        std::vector<uint8_t> values_uint8_t(values.begin(), values.end());
        return cube(variables, values_uint8_t);
    }

    SylvanUBDD cube(std::vector<SylvanUBDD> variables) const override {
        BDDsylvan bdd = bdd_.VectorCube(extractBDDs(variables));
        return SylvanUBDD(bdd);
    }

    SylvanUBDD permute(const std::vector<size_t> &from, const std::vector<size_t> &to) const override {
        return SylvanUBDD(bdd_.Permute(std::vector<uint32_t>(from.begin(), from.end()), std::vector<uint32_t>(to.begin(), to.end())));
    }

    double countMinterm(size_t nvars) const override {
        size_t nodeCount = bdd_.NodeCount();
        if (nodeCount <= nvars) {
            return bdd_.SatCount(nvars);
        } else {
            return bdd_.SatCount(nodeCount) / std::pow(2.0, nodeCount - nvars);
        }
    }

    void printMinterm() const override {
        sylvan_print(bdd_.GetBDD());
    }


    SylvanUBDD existAbstract(const SylvanUBDD &cube) const override {
        return SylvanUBDD(bdd_.ExistAbstract(cube.bdd_));
    }

    SylvanUBDD andAbstract(const SylvanUBDD &g, const SylvanUBDD &cube) const override {
        BDDsylvan bdd = bdd_.AndAbstract(g.bdd_, cube.bdd_);
        return SylvanUBDD(bdd);
    }

    UBDDMintermIterator *generateMintermIterator(std::vector<size_t> &ivars) const override {
        BDDsylvan bddSet = var(ivars[0]).bdd_;
        for (size_t i = 1; i < ivars.size(); i++) {
            bddSet &= var(ivars[i]).bdd_;
        }
        return new SylvanUBDDMintermIterator(bdd_.GetBDD(), bddSet.GetBDD(), ivars);
    }

    SylvanUBDD complement(SylvanUBDD &symbolicSet,
                          std::vector<size_t> nofGridPoints,
                          std::vector<size_t> nofBddVars,
                          std::vector<std::vector<size_t>> indBddVars,
                          size_t dim) override {
        /* generate ADD variables and ADD for the grid points in each dimension */
        Mtbdd constant;
        std::vector<Mtbdd> aVar(dim);
        std::vector<std::vector<Mtbdd>> addVars(dim);
        for (size_t i = 0; i < dim; i++) {
            aVar[i] = Mtbdd::mtbddZero();
            addVars[i] = std::vector<Mtbdd>(nofBddVars[i]);
            double c = 1;
            for (size_t j = 0; j < nofBddVars[i]; j++) {
                constant = Mtbdd::doubleTerminal(c);
                addVars[i][j] = Mtbdd::mtbddVar(indBddVars[i][j]); /* the argument should be uint32_t */
                addVars[i][j] *= constant;
                aVar[i] += addVars[i][j];
                c *= 2;
            }
        }
        /* set grid points outside of the uniform grid to infinity */
        BDDsylvan outside = zero().bdd_;
        for (size_t i = 0; i < dim; i++) {
            BDDsylvan bdd = aVar[i].BddThreshold(nofGridPoints[i]);
            outside += bdd;
        }

        BDDsylvan inside = !outside;
        return SylvanUBDD(inside) & !symbolicSet;
    }

    /* compute symbolic representation (=bdd) of a polytope { x | Hx<= h }*/
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
                               const std::vector<size_t> &nofGridPoints) override {
        /* define minusInfinity (for outside area of the domain) and epsilon */
        double unusedFloat = firstGridPoint[0] - eta[0];
        for (size_t i = 1; i < dim; i++) {
            if (unusedFloat > firstGridPoint[i] - eta[i]) {
                unusedFloat = firstGridPoint[i] - eta[i];
            }
        }
        double minusInfinity = unusedFloat;
        double epsilon = DBL_MIN;
        /* generate ADD variables and ADD for the grid points in each dimension */
        const int OUTER = 0;
        const int INNER = 1;
        Mtbdd constant;
        std::vector<Mtbdd> aVar(dim);
        std::vector<std::vector<Mtbdd>> addVars(dim);
        for (size_t i = 0; i < dim; i++) {
            aVar[i] = Mtbdd::mtbddZero();
            addVars[i] = std::vector<Mtbdd>(nofBddVars[i]);
            double c = 1.0;
            for (size_t j = 0; j < nofBddVars[i]; j++) {
                constant = Mtbdd::doubleTerminal(c);
                addVars[i][j] = Mtbdd::mtbddVar(indBddVars[i][j]); /* the argument should be uint32_t */
                addVars[i][j] *= constant;
                aVar[i] += addVars[i][j];
                c *= 2;
            }
        }
        /* set grid points outside of the uniform grid to infinity */
        auto outside = Mtbdd::mtbddZero();
        for (size_t i = 0; i < dim; i++) {
            Mtbdd bdd = aVar[i].BddThreshold(nofGridPoints[i]);
            Mtbdd add = bdd * minusInfinity;
            outside += add;
        }

        /* compute values for polytope */
        BDDsylvan polytope = Bdd::bddOne();
        for (size_t i = 0; i < p; i++) {
            /* update set of outside nodes */
            Mtbdd bdd = !polytope;
            outside += bdd * minusInfinity;
            Mtbdd halfspace = outside;
            for (size_t j = 0; j < dim; j++) {
                if (H[i * dim + j] == 0)
                    continue;
                Mtbdd c = Mtbdd::doubleTerminal((-H[i * dim + j]));
                Mtbdd addEta = Mtbdd::doubleTerminal(eta[j]);
                Mtbdd radius = Mtbdd::doubleTerminal((eta[j] / 2.0 + z[j]));
                Mtbdd first = Mtbdd::doubleTerminal(firstGridPoint[j]);
                if ((H[i * dim + j] >= 0 && type == OUTER) || (H[i * dim + j] < 0 && type == INNER))
                    halfspace += ((aVar[j] * addEta) + (first - radius)) * c;
                else
                    halfspace += ((aVar[j] * addEta) + (first + radius)) * c;
            }
            /* add halfspace to polytope */
            if (type == OUTER)
                polytope &= halfspace.BddThreshold(-(h[i] + epsilon));
            else
                polytope &= halfspace.BddThreshold(-(h[i] - epsilon));
        }
        return SylvanUBDD(polytope);
    }

    int save(FILE *file) override {
        MTBDD m = bdd_.GetBDD();
        mtbdd_writer_tobinary(file, &m, 1);
        return 1;
    }

    SylvanUBDD load(FILE *file, std::vector<int> composeids, int newID) override {
        MTBDD m;
        mtbdd_reader_frombinary(file, &m, 1);
        BDDsylvan bdd(m);
        SylvanUBDD ubdd(bdd);
        if (newID) {
            std::vector<size_t> from(composeids.size());
            std::vector<size_t> to(composeids.begin(), composeids.end());
            for (int i = 0; i < composeids.size(); i++)
                from[i] = i;
            return ubdd.permute(from, to);
        }
        return ubdd;
    }

    SylvanUBDD transfer(const SylvanUBDD &destination) const override {
        return *this;
    }

    bool isCoverEqual(const SylvanUBDD &other) const override { return true; }

    SylvanUBDD &operator=(const SylvanUBDD &right) {
        bdd_ = right.bdd_;
        return *this;
    }
    bool operator==(const SylvanUBDD &other) const override { return bdd_ == other.bdd_; }
    bool operator!=(const SylvanUBDD &other) const override { return bdd_ != other.bdd_; }
    bool operator<=(const SylvanUBDD &other) const override { return bdd_ <= other.bdd_; }
    bool operator>=(const SylvanUBDD &other) const override { return bdd_ >= other.bdd_; }
    bool operator<(const SylvanUBDD &other) const override { return bdd_ < other.bdd_; }
    bool operator>(const SylvanUBDD &other) const override { return bdd_ > other.bdd_; }
    SylvanUBDD operator!() const override {
        BDDsylvan bdd = !bdd_;
        return SylvanUBDD(bdd);
    }
    SylvanUBDD operator~() const override {
        BDDsylvan bdd = ~bdd_;
        return SylvanUBDD(bdd);
    }
    SylvanUBDD operator*(const SylvanUBDD &other) const override {
        BDDsylvan bdd = bdd_ * other.bdd_;
        return SylvanUBDD(bdd);
    }
    SylvanUBDD operator*=(const SylvanUBDD &other) override {
        bdd_ *= other.bdd_;
        return *this;
    }
    SylvanUBDD operator&(const SylvanUBDD &other) const override {
        BDDsylvan bdd = bdd_ & other.bdd_;
        return SylvanUBDD(bdd);
    }
    SylvanUBDD operator&=(const SylvanUBDD &other) override {
        bdd_ &= other.bdd_;
        return *this;
    }
    SylvanUBDD operator+(const SylvanUBDD &other) const override {
        BDDsylvan bdd = bdd_ + other.bdd_;
        return SylvanUBDD(bdd);
    }
    SylvanUBDD operator+=(const SylvanUBDD &other) override {
        bdd_ += other.bdd_;
        return *this;
    }
    SylvanUBDD operator|(const SylvanUBDD &other) const override {
        BDDsylvan bdd = bdd_ | other.bdd_;
        return SylvanUBDD(bdd);
    }
    SylvanUBDD operator|=(const SylvanUBDD &other) override {
        bdd_ |= other.bdd_;
        return *this;
    }
    SylvanUBDD operator^(const SylvanUBDD &other) const override {
        BDDsylvan bdd = bdd_ ^ other.bdd_;
        return SylvanUBDD(bdd);
    }
    SylvanUBDD operator^=(const SylvanUBDD &other) override {
        bdd_ ^= other.bdd_;
        return *this;
    }
    SylvanUBDD operator-(const SylvanUBDD &other) const override {
        BDDsylvan bdd = bdd_ - other.bdd_;
        return SylvanUBDD(bdd);
    }
    SylvanUBDD operator-=(const SylvanUBDD &other) override {
        bdd_ -= other.bdd_;
        return *this;
    }

private:
    static std::vector<BDDsylvan> extractBDDs(const std::vector<SylvanUBDD> &ubdds) {
        std::vector<BDDsylvan> bdds(ubdds.size());
        for (size_t i = 0; i < ubdds.size(); i++)
            bdds[i] = ubdds[i].bdd_;
        return bdds;
    }
};
//std::map<size_t, SylvanUBDD *> SylvanUBDD::nodes_map = SylvanUBDD::create_map();

#endif /* SYLVANUBDD_HH_
 */
