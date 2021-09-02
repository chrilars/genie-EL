/*
 * CuddUBDD.hh
 *
 *  created on: 01.07.2021 (In Progress)
 *      author: Mateusz Rychlicki
 *
 */
#include "BaseUBDD.hh"
#include "CuddUBDDMintermIterator.hh"
#include <complex>
#include <memory>
#include <utility>

#include "cudd.h"
#include "cuddObj.hh"
#include "dddmp.h"

#ifndef CUDDUBDD_HH_
#    define CUDDUBDD_HH_

typedef BDD BDDcudd;

class CuddUBDD : public BaseUBDD<CuddUBDD> {
public:
    BDDcudd bdd_;
    std::shared_ptr<Cudd> cudd_;

    CuddUBDD(BDDcudd &bdd, std::shared_ptr<Cudd> cudd) {
        cudd_ = cudd;
        bdd_ = bdd;
    }

    CuddUBDD(unsigned int numVars,
             unsigned int numVarsZ,
             unsigned int numSlots = CUDD_UNIQUE_SLOTS,
             unsigned int cacheSize = CUDD_CACHE_SLOTS,
             unsigned long maxMemory = 0,
             PFC defaultHandler = defaultError) {
        cudd_ = std::make_shared<Cudd>(numVars, numVarsZ, numSlots, cacheSize, maxMemory, defaultHandler);
    }

    CuddUBDD(const Cudd &x) {
        cudd_ = std::make_shared<Cudd>(x);
    }

    CuddUBDD(CuddUBDD &ubdd) {
        cudd_ = ubdd.cudd_;
        bdd_ = BDD(ubdd.bdd_);
    }

    CuddUBDD(const CuddUBDD &ubdd) {
        cudd_ = ubdd.cudd_;
        bdd_ = BDD(ubdd.bdd_);
    }

    CuddUBDD() {
        cudd_ = std::make_shared<Cudd>();
    }

    CuddUBDD zero() const override {
        BDDcudd bdd = cudd_->bddZero();
        return CuddUBDD(bdd, cudd_);
    }

    CuddUBDD one() const override {
        BDDcudd bdd = cudd_->bddOne();
        return CuddUBDD(bdd, cudd_);
    }

    CuddUBDD var() const override {
        BDDcudd bdd = cudd_->bddVar();
        return CuddUBDD(bdd, cudd_);
    }

    CuddUBDD var(size_t index) const override {
        BDDcudd bdd = cudd_->bddVar(index);
        return CuddUBDD(bdd, cudd_);
    }

    size_t nodeSize() const override {
        return cudd_->ReadSize();
    }

    size_t nodeIndex() const override {
        return bdd_.NodeReadIndex();
    }

    size_t nodeCount() const {
        return bdd_.SupportSize();
    }

    CuddUBDD cube(std::vector<CuddUBDD> variables, std::vector<int> values) const override {
        std::vector<BDDcudd> varsBdd = extractBDDs(variables);
        BDDcudd bdd = variables[0].cudd_->bddComputeCube(&varsBdd[0], &values[0], variables.size());
        return CuddUBDD(bdd, variables[0].cudd_);
    }

    CuddUBDD cube(std::vector<CuddUBDD> variables, std::vector<uint8_t> values) const override {
        std::vector<int> phase_int(values.begin(), values.end());
        return cube(variables, phase_int);
    }

    CuddUBDD cube(std::vector<CuddUBDD> variables) const override {
        std::vector<BDDcudd> varsBdd = extractBDDs(variables);
        BDDcudd bdd = cudd_->bddComputeCube(&varsBdd[0], nullptr, variables.size());// Same as BDD::computeCube
        return CuddUBDD(bdd, cudd_);
    }

    CuddUBDD permute(const std::vector<size_t> &from, const std::vector<size_t> &to) const override {
        BDDcudd bdd;
        std::vector<int> permutation;
        if (nodeSize() > from.size()) {
            permutation = std::vector<int>(nodeSize());
            for (size_t i = 0; i < permutation.size(); i++)
                permutation[i] = i;
        } else {
            permutation = std::vector<int>(from.size());
        }
        for (size_t i = 0; i < from.size(); i++)
            permutation[from[i]] = to[i];

        bdd = bdd_.Permute(&permutation[0]);
        return CuddUBDD(bdd, cudd_);
    }

    double countMinterm(size_t nvars) const override {
        return bdd_.CountMinterm(nvars);
    }

    void printMinterm() const override {
        bdd_.PrintMinterm();
    }

    double readEpsilon() const {
        return cudd_->ReadEpsilon();
    }


    CuddUBDD existAbstract(const CuddUBDD &cube) const override {
        BDD bdd = bdd_.ExistAbstract(cube.bdd_);
        return CuddUBDD(bdd, cudd_);
    }

    CuddUBDD andAbstract(const CuddUBDD &g, const CuddUBDD &cube) const override {
        BDDcudd bdd = bdd_.AndAbstract(g.bdd_, cube.bdd_);
        return CuddUBDD(bdd, g.cudd_);
    }

    UBDDMintermIterator *generateMintermIterator(std::vector<size_t> &ivars) const override {
        return new CuddUBDDMintermIterator(bdd_, ivars);
    }

    CuddUBDD complement(CuddUBDD &symbolicSet,
                        std::vector<size_t> nofGridPoints,
                        std::vector<size_t> nofBddVars,
                        std::vector<std::vector<size_t>> indBddVars,
                        size_t dim) override {
        /* generate ADD variables and ADD for the grid points in each dimension */
        ADD constant;
        std::vector<ADD> aVar(dim);
        std::vector<std::vector<ADD>> addVars(dim);
        for (size_t i = 0; i < dim; i++) {
            aVar[i] = cudd_->addZero();
            addVars[i] = std::vector<ADD>(nofBddVars[i]);
            CUDD_VALUE_TYPE c = 1;
            for (size_t j = 0; j < nofBddVars[i]; j++) {
                constant = cudd_->constant(c);
                addVars[i][j] = cudd_->addVar(indBddVars[i][j]);
                addVars[i][j] *= constant;
                aVar[i] += addVars[i][j];
                c *= 2;
            }
        }
        /* set grid points outside of the uniform grid to infinity */
        BDDcudd outside = cudd_->bddZero();
        for (size_t i = 0; i < dim; i++) {
            BDDcudd bdd = aVar[i].BddThreshold(nofGridPoints[i]);
            outside += bdd;
        }

        BDDcudd inside = !outside;
        return CuddUBDD(inside, cudd_) & !symbolicSet;
    }

    /* compute symbolic representation (=bdd) of a polytope { x | Hx<= h }*/
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
                             const std::vector<size_t> &nofGridPoints) override {
        /* generate ADD variables and ADD for the grid points in each dimension */
        const int OUTER = 0;
        const int INNER = 1;
        ADD constant;
        ADD *aVar = new ADD[dim];
        ADD **addVars = new ADD *[dim];
        for (size_t i = 0; i < dim; i++) {
            aVar[i] = cudd_->addZero();
            addVars[i] = new ADD[nofBddVars[i]];
            CUDD_VALUE_TYPE c = 1;
            for (size_t j = 0; j < nofBddVars[i]; j++) {
                constant = cudd_->constant(c);
                addVars[i][j] = cudd_->addVar(indBddVars[i][j]);
                addVars[i][j] *= constant;
                aVar[i] += addVars[i][j];
                c *= 2;
            }
        }
        /* set grid points outside of the uniform grid to infinity */
        ADD outside = cudd_->addZero();
        for (size_t i = 0; i < dim; i++) {
            BDDcudd bdd = aVar[i].BddThreshold(nofGridPoints[i]);
            ADD add = bdd.Add() * cudd_->minusInfinity();
            outside += add;
        }

        /* compute values for polytope */
        BDDcudd polytope = cudd_->bddOne();
        for (size_t i = 0; i < p; i++) {
            /* update set of outside nodes */
            BDDcudd bdd = !polytope;
            outside += bdd.Add() * cudd_->minusInfinity();
            ADD halfspace = outside;
            for (size_t j = 0; j < dim; j++) {
                if (H[i * dim + j] == 0)
                    continue;
                ADD c = cudd_->constant((-H[i * dim + j]));
                ADD addEta = cudd_->constant(eta[j]);
                ADD radius = cudd_->constant((eta[j] / 2.0 + z[j]));
                ADD first = cudd_->constant(firstGridPoint[j]);
                if ((H[i * dim + j] >= 0 && type == OUTER) || (H[i * dim + j] < 0 && type == INNER))
                    halfspace += ((aVar[j] * addEta) + (first - radius)) * c;
                else
                    halfspace += ((aVar[j] * addEta) + (first + radius)) * c;
            }
            /* add halfspace to polytope */
            if (type == OUTER)
                polytope &= halfspace.BddThreshold(-(h[i] + cudd_->ReadEpsilon()));
            else
                polytope &= halfspace.BddThreshold(-(h[i] - cudd_->ReadEpsilon()));
        }
        delete[] aVar;
        for (size_t i = 0; i < dim; i++)
            delete[] addVars[i];
        delete[] addVars;
        return CuddUBDD(polytope, cudd_);
    }

    int save(FILE *file) override{
        /* before we save the BDD to file, we save it to another manager,
         * because the current manager is classified as ADD manager */
        Cudd mdest;
        BDD tosave = bdd_.Transfer(mdest);
        return Dddmp_cuddBddStore(
                mdest.getManager(),
                nullptr,
                tosave.getNode(),
                nullptr,
                nullptr,
                DDDMP_MODE_BINARY,
                DDDMP_VARIDS,
                nullptr,
                file);
    }

    CuddUBDD load(FILE *file, std::vector<int> composeids, int newID) override{
        /* if match then we have to create new variable id's and load the bdd with those new ids */
        DdNode *dd = Dddmp_cuddBddLoad(cudd_->getManager(),
                                       (newID) ? DDDMP_VAR_COMPOSEIDS : DDDMP_VAR_MATCHIDS,
                                       nullptr,
                                       nullptr,
                                       &composeids[0],
                                       DDDMP_MODE_BINARY,
                                       nullptr,
                                       file);
        BDDcudd bdd(*cudd_, dd);
        return CuddUBDD(bdd, cudd_);
    }

    CuddUBDD transfer(const CuddUBDD &destination) const override{
        BDD bdd = BDD(bdd_);
        return CuddUBDD(bdd, cudd_);
    }

    bool isCoverEqual(const CuddUBDD &other) const override{
        return cudd_->getManager() == other.cudd_->getManager();
    }

    CuddUBDD& operator=(const CuddUBDD &right) override {
        bdd_ = right.bdd_;
        cudd_ = right.cudd_;
        return *this;
    }
    bool operator==(const CuddUBDD &other) const override { return bdd_ == other.bdd_; }
    bool operator!=(const CuddUBDD &other) const override { return bdd_ != other.bdd_; }
    bool operator<=(const CuddUBDD &other) const override { return bdd_ <= other.bdd_; }
    bool operator>=(const CuddUBDD &other) const override { return bdd_ >= other.bdd_; }
    bool operator<(const CuddUBDD &other) const override { return bdd_ < other.bdd_; }
    bool operator>(const CuddUBDD &other) const override { return bdd_ > other.bdd_; }
    CuddUBDD operator!() const override {
        BDDcudd bdd = !bdd_;
        return CuddUBDD(bdd, cudd_);
    }
    CuddUBDD operator~() const override {
        BDDcudd bdd = ~bdd_;
        return CuddUBDD(bdd, cudd_);
    }
    CuddUBDD operator*(const CuddUBDD &other) const override {
        BDDcudd bdd = bdd_ * other.bdd_;
        return CuddUBDD(bdd, cudd_);
    }
    CuddUBDD operator*=(const CuddUBDD &other) override {
        bdd_ *= other.bdd_;
        return *this;
    }
    CuddUBDD operator&(const CuddUBDD &other) const override {
        BDDcudd bdd = bdd_ & other.bdd_;
        return CuddUBDD(bdd, cudd_);
    }
    CuddUBDD operator&=(const CuddUBDD &other) override {
        bdd_ &= other.bdd_;
        return *this;
    }
    CuddUBDD operator+(const CuddUBDD &other) const override {
        BDDcudd bdd = bdd_ + other.bdd_;
        return CuddUBDD(bdd, cudd_);
    }
    CuddUBDD operator+=(const CuddUBDD &other) override {
        bdd_ += other.bdd_;
        return *this;
    }
    CuddUBDD operator|(const CuddUBDD &other) const override {
        BDDcudd bdd = bdd_ | other.bdd_;
        return CuddUBDD(bdd, cudd_);
    }
    CuddUBDD operator|=(const CuddUBDD &other) override {
        bdd_ |= other.bdd_;
        return *this;
    }
    CuddUBDD operator^(const CuddUBDD &other) const override {
        BDDcudd bdd = bdd_ ^ other.bdd_;
        return CuddUBDD(bdd, cudd_);
    }
    CuddUBDD operator^=(const CuddUBDD &other) override {
        bdd_ ^= other.bdd_;
        return *this;
    }
    CuddUBDD operator-(const CuddUBDD &other) const override {
        BDDcudd bdd = bdd_ - other.bdd_;
        return CuddUBDD(bdd, cudd_);
    }
    CuddUBDD operator-=(const CuddUBDD &other) override {
        bdd_ -= other.bdd_;
        return *this;
    }

    /* compute the smallest radius over all 2-norm balls that contain
     * the set  L([-eta[0]/2, eta[0]2]x ... x [eta[dim-1]/2, eta[dim-1]/2]) */
    double computeEllipsoidRadius(const std::vector<double> L, size_t dim, const std::vector<double> eta, const std::vector<double> z) const {

        /* compute smallest radius of ball that contains L*cell */
        Cudd mgr;
        ADD *y = new ADD[dim];
        ADD *x = new ADD[dim];
        ADD *v = new ADD[dim];
        for (size_t i = 0; i < dim; i++) {
            v[i] = mgr.addVar(i);
            x[i] = v[i] * mgr.constant(eta[i]) - mgr.constant(eta[i] / 2.0 + z[i]);
        }
        /* compute y= Lx */
        ADD lij;
        for (size_t i = 0; i < dim; i++) {
            y[i] = mgr.addZero();
            for (size_t j = 0; j < dim; j++) {
                lij = mgr.constant(L[dim * i + j]);
                y[i] += lij * x[j];
            }
        }
        /* compute r = x'L'Lx */
        ADD r = mgr.addZero();
        for (size_t i = 0; i < dim; i++)
            r += y[i] * y[i];
        ADD max = r.FindMax();
        delete[] x;
        delete[] y;
        delete[] v;

        return std::sqrt(Cudd_V(max.getNode())) + cudd_->ReadEpsilon();
    }

private:
    static std::vector<BDDcudd> extractBDDs(std::vector<CuddUBDD> ubdds) {
        std::vector<BDDcudd> bdds(ubdds.size());
        for (size_t i = 0; i < ubdds.size(); i++)
            bdds[i] = ubdds[i].bdd_;
        return bdds;
    }
};

#endif /* CUDDUBDD_HH_
 */
