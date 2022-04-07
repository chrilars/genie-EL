/** @file CuddUBDD.cpp
*
*  @date 10.09.2021
*  @author Mateusz Rychlicki
*/

#include "ubdd/CuddUBDD.hh"
#include "ubdd/CuddUBDDMintermIterator.hh"
#include "ubdd/BaseUBDD.hh"
#include <complex>
#include <memory>
#include <utility>
#include <vector>
//duap
#include "cudd.h"
#include "cuddObj.hh"
#include "dddmp.h"

namespace fairsyn {

    CuddUBDD::CuddUBDD(BDDcudd &bdd, std::shared_ptr<Cudd> cudd) {
        cudd_ = cudd;
        bdd_ = bdd;
    }

    CuddUBDD::CuddUBDD(unsigned int numVars,
                       unsigned int numVarsZ,
                       unsigned int numSlots,
                       unsigned int cacheSize,
                       unsigned long maxMemory,
                       PFC defaultHandler) {
        cudd_ = std::make_shared<Cudd>(numVars, numVarsZ, numSlots, cacheSize, maxMemory, defaultHandler);
    }

    CuddUBDD::CuddUBDD(const Cudd &x) {
        cudd_ = std::make_shared<Cudd>(x);
    }

    CuddUBDD::CuddUBDD(CuddUBDD &ubdd) {
        cudd_ = ubdd.cudd_;
        bdd_ = BDD(ubdd.bdd_);
    }

    CuddUBDD::CuddUBDD(const CuddUBDD &ubdd) {
        cudd_ = ubdd.cudd_;
        bdd_ = BDD(ubdd.bdd_);
    }

    CuddUBDD::CuddUBDD() {
        cudd_ = std::make_shared<Cudd>();
    }

    CuddUBDD CuddUBDD::zero() const {
        BDDcudd bdd = cudd_->bddZero();
        return CuddUBDD(bdd, cudd_);
    }

    CuddUBDD CuddUBDD::one() const {
        BDDcudd bdd = cudd_->bddOne();
        return CuddUBDD(bdd, cudd_);
    }

    CuddUBDD CuddUBDD::var() const {
        BDDcudd bdd = cudd_->bddVar();
        return CuddUBDD(bdd, cudd_);
    }

    CuddUBDD CuddUBDD::var(size_t index) const {
        BDDcudd bdd = cudd_->bddVar(index);
        return CuddUBDD(bdd, cudd_);
    }

    size_t CuddUBDD::nodeSize() const {
        return cudd_->ReadSize();
    }

    size_t CuddUBDD::nodeIndex() const {
        return bdd_.NodeReadIndex();
    }

    size_t CuddUBDD::nodeCount() const {
        return bdd_.SupportSize();
    }

    CuddUBDD CuddUBDD::cube(std::vector<CuddUBDD> variables, std::vector<int> values) const {
        std::vector<BDDcudd> varsBdd = extractBDDs(variables);
        BDDcudd bdd = variables[0].cudd_->bddComputeCube(&varsBdd[0], &values[0], variables.size());
        return CuddUBDD(bdd, variables[0].cudd_);
    }

    CuddUBDD CuddUBDD::cube(std::vector<CuddUBDD> variables, std::vector<uint8_t> values) const {
        std::vector<int> phase_int(values.begin(), values.end());
        return cube(variables, phase_int);
    }

    CuddUBDD CuddUBDD::cube(std::vector<CuddUBDD> variables) const {
        std::vector<BDDcudd> varsBdd = extractBDDs(variables);
        BDDcudd bdd = cudd_->bddComputeCube(&varsBdd[0], nullptr, variables.size());// Same as BDD::computeCube
        return CuddUBDD(bdd, cudd_);
    }

    CuddUBDD CuddUBDD::permute(const std::vector<size_t> &from, const std::vector<size_t> &to) const {
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

    double CuddUBDD::countMinterm(size_t nvars) const {
        return bdd_.CountMinterm(nvars);
    }

    void CuddUBDD::printMinterm() const {
        bdd_.PrintMinterm();
    }

    double CuddUBDD::readEpsilon() const {
        return cudd_->ReadEpsilon();
    }


    CuddUBDD CuddUBDD::existAbstract(const CuddUBDD &cube) const {
        BDD bdd = bdd_.ExistAbstract(cube.bdd_);
        return CuddUBDD(bdd, cudd_);
    }

    CuddUBDD CuddUBDD::andAbstract(const CuddUBDD &g, const CuddUBDD &cube) const {
        BDDcudd bdd = bdd_.AndAbstract(g.bdd_, cube.bdd_);
        return CuddUBDD(bdd, g.cudd_);
    }

    UBDDMintermIterator *CuddUBDD::generateMintermIterator(std::vector<size_t> &ivars) const {
        return new CuddUBDDMintermIterator(bdd_, ivars);
    }

    CuddUBDD CuddUBDD::complement(CuddUBDD &symbolicSet,
                                  std::vector<size_t> nofGridPoints,
                                  std::vector<size_t> nofBddVars,
                                  std::vector<std::vector<size_t>> indBddVars,
                                  size_t dim) {
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
    CuddUBDD CuddUBDD::computePolytope(const size_t p,
                                       const std::vector<double> &H,
                                       const std::vector<double> &h,
                                       int type,
                                       size_t dim,
                                       const std::vector<double> &eta,
                                       const std::vector<double> &z,
                                       const std::vector<double> &firstGridPoint,
                                       const std::vector<size_t> &nofBddVars,
                                       const std::vector<std::vector<size_t>> &indBddVars,
                                       const std::vector<size_t> &nofGridPoints) {
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

    int CuddUBDD::save(FILE *file) {
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

    CuddUBDD CuddUBDD::load(FILE *file, std::vector<int> composeids, int newID) {
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

    CuddUBDD CuddUBDD::transfer(const CuddUBDD &destination) const {
        BDD bdd = BDD(bdd_);
        return CuddUBDD(bdd, cudd_);
    }

    bool CuddUBDD::isCoverEqual(const CuddUBDD &other) const {
        return cudd_->getManager() == other.cudd_->getManager();
    }

    CuddUBDD &CuddUBDD::operator=(const CuddUBDD &right) {
        bdd_ = right.bdd_;
        cudd_ = right.cudd_;
        return *this;
    }
    bool CuddUBDD::operator==(const CuddUBDD &other) const { return bdd_ == other.bdd_; }
    bool CuddUBDD::operator!=(const CuddUBDD &other) const { return bdd_ != other.bdd_; }
    bool CuddUBDD::operator<=(const CuddUBDD &other) const { return bdd_ <= other.bdd_; }
    bool CuddUBDD::operator>=(const CuddUBDD &other) const { return bdd_ >= other.bdd_; }
    bool CuddUBDD::operator<(const CuddUBDD &other) const { return bdd_ < other.bdd_; }
    bool CuddUBDD::operator>(const CuddUBDD &other) const { return bdd_ > other.bdd_; }
    CuddUBDD CuddUBDD::operator!() const {
        BDDcudd bdd = !bdd_;
        return CuddUBDD(bdd, cudd_);
    }
    CuddUBDD CuddUBDD::operator~() const {
        BDDcudd bdd = ~bdd_;
        return CuddUBDD(bdd, cudd_);
    }
    CuddUBDD CuddUBDD::operator*(const CuddUBDD &other) const {
        BDDcudd bdd = bdd_ * other.bdd_;
        return CuddUBDD(bdd, cudd_);
    }
    CuddUBDD CuddUBDD::operator*=(const CuddUBDD &other) {
        bdd_ *= other.bdd_;
        return *this;
    }
    CuddUBDD CuddUBDD::operator&(const CuddUBDD &other) const {
        BDDcudd bdd = bdd_ & other.bdd_;
        return CuddUBDD(bdd, cudd_);
    }
    CuddUBDD CuddUBDD::operator&=(const CuddUBDD &other) {
        bdd_ &= other.bdd_;
        return *this;
    }
    CuddUBDD CuddUBDD::operator+(const CuddUBDD &other) const {
        BDDcudd bdd = bdd_ + other.bdd_;
        return CuddUBDD(bdd, cudd_);
    }
    CuddUBDD CuddUBDD::operator+=(const CuddUBDD &other) {
        bdd_ += other.bdd_;
        return *this;
    }
    CuddUBDD CuddUBDD::operator|(const CuddUBDD &other) const {
        BDDcudd bdd = bdd_ | other.bdd_;
        return CuddUBDD(bdd, cudd_);
    }
    CuddUBDD CuddUBDD::operator|=(const CuddUBDD &other) {
        bdd_ |= other.bdd_;
        return *this;
    }
    CuddUBDD CuddUBDD::operator^(const CuddUBDD &other) const {
        BDDcudd bdd = bdd_ ^ other.bdd_;
        return CuddUBDD(bdd, cudd_);
    }
    CuddUBDD CuddUBDD::operator^=(const CuddUBDD &other) {
        bdd_ ^= other.bdd_;
        return *this;
    }
    CuddUBDD CuddUBDD::operator-(const CuddUBDD &other) const {
        BDDcudd bdd = bdd_ - other.bdd_;
        return CuddUBDD(bdd, cudd_);
    }
    CuddUBDD CuddUBDD::operator-=(const CuddUBDD &other) {
        bdd_ -= other.bdd_;
        return *this;
    }

    /* compute the smallest radius over all 2-norm balls that contain
     * the set  L([-eta[0]/2, eta[0]2]x ... x [eta[dim-1]/2, eta[dim-1]/2]) */
    double CuddUBDD::computeEllipsoidRadius(const std::vector<double> L, size_t dim, const std::vector<double> eta, const std::vector<double> z) const {

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
}// namespace fairsyn