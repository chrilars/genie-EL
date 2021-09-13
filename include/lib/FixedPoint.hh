/*
* FixedPoint.hh
*
*/

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>
#include <sylvan_obj.hpp>
#include <vector>

#include "lib/Arena.hh"
#include "lib/RabinAutomaton.hh"
#include "utils/TicToc.hh"

struct const_arg_recursive_rabin {
    const bool accl_on;
    const size_t M; /* the bound on the iteration count for memorizing the BDDs from the past iterations */
    const int depth;
    const std::vector<fairsyn::rabin_pair_> pairs;
    const sylvan::Bdd initial_seed;
    const int verbose;
};

struct nconst_arg_recursive_rabin {
    sylvan::Bdd seqR;
    sylvan::Bdd right;
    std::vector<size_t> *indexRP;
    std::vector<size_t> *indexY;
    std::vector<size_t> *indexX;

    //////////////////// for acceleration
    std::map<std::vector<size_t>, size_t> *permutation;
    std::vector<std::map<std::vector<size_t>, sylvan::Bdd>> *hist_Y;
    std::vector<std::map<std::vector<size_t>, sylvan::Bdd>> *hist_X;
};

namespace fairsyn {
    /*
    * class: FixedPoint
    *
    *
    * provides the fixed point computation for the Rabin specification
    * on finite transition systems with the edge fairness condition
    *
    */
    class FixedPoint : public Arena {
    public:
        /* var: tr_domain_
    *transitoin bdd with cubePost_ abstracted */
        sylvan::Bdd tr_domain_;
        /* var: live_domain
    * set of live edges with cubePost_ abstracted */
        sylvan::Bdd live_domain_;
        /* cubes with post variables; used in the existential abstraction  */
        sylvan::Bdd cubePost_;
        /* cubes with other variables (outside the pre and the post variables) on which the transitions possibly depend */
        sylvan::Bdd cubeOther_;
        /* var: RabinPairs_
    * vector of the rabin pairs*/
        std::vector<rabin_pair_> RabinPairs_;

    public:
        /* constructor: FixedPoint
    *
    * The rabin pairs are given as a vector of the pairs (Gi,Ri) for the specification (GF Gi & FG !Ri)
    *
    */
        FixedPoint(std::vector<size_t> &nodes,
                   std::vector<size_t> &sys_nodes,
                   std::vector<size_t> &env_nodes,
                   std::vector<std::vector<size_t>> &transitions,
                   std::vector<std::vector<size_t>> &live_edges,
                   std::vector<std::array<std::vector<size_t>, 2>> &RabinPairs, std::set<uint32_t> &all_variables) : Arena(all_variables, nodes, sys_nodes, env_nodes, transitions, live_edges) {
            cubePost_ = sylvan::Bdd::VariablesCube(sort(postVars_));
            cubeOther_ = sylvan::Bdd::bddOne();
            tr_domain_ = tr_.ExistAbstract(cubePost_);
            live_domain_ = live_.ExistAbstract(cubePost_);
            sylvan::BddSet preBddVars = sylvan::BddSet::fromVector(preVars_);
            for (size_t i = 0; i < RabinPairs.size(); ++i) {
                rabin_pair_ p;
                p.rabin_index_ = i;
                p.G_ = SetToBdd(RabinPairs[i][0], preBddVars, preVars_.size());
                p.nR_ = (!(SetToBdd(RabinPairs[i][1], preBddVars, preVars_.size()))) & nodes_;
                RabinPairs_.push_back(p);
            }
        }
        /*
       * constructor
       */
        FixedPoint(const Arena *other,
                   const RabinAutomaton *rabin) {
            /* the set of nodes of this automaton is the product (computed using BDD product) state space */
            nodes_ = other->nodes_ * rabin->stateSpace_;
            /* the system nodes in the product are those nodes whose respective component in "other" belongs to the system nodes in "other" */
            sys_nodes_ = other->sys_nodes_;
            /* similarly, the environment nodes */
            env_nodes_ = other->env_nodes_;
            /* the set of variables is the union of variables of the arena and the rabin automaton */
            preVars_ = other->preVars_;
            for (auto i = rabin->stateVars_.begin(); i != rabin->stateVars_.end(); ++i) {
                preVars_.push_back(*i);
            }
            postVars_ = other->postVars_;
            for (auto i = rabin->postStateVars_.begin(); i != rabin->postStateVars_.end(); ++i) {
                postVars_.push_back(*i);
            }
            /* the joint transitions are the product of the two respective transitions */
            tr_ = other->tr_ * rabin->transitions_;
            live_ = other->live_ * rabin->transitions_;
            /* compute the helper BDDs for synthesis */
            cubePost_ = sylvan::Bdd::VariablesCube(sort(postVars_));
            cubeOther_ = sylvan::Bdd::bddOne();
            tr_domain_ = tr_.ExistAbstract(cubePost_);
            live_domain_ = live_.ExistAbstract(cubePost_);
            /* the rabin pairs are directly inherited from the rabin automaton */
            RabinPairs_ = rabin->RabinPairs_;
        }
        /*
       * constructor
       */
        FixedPoint(const Arena *other,
                   const RabinAutomaton *rabin,
                   const std::vector<uint32_t> other_variables) {
            /* the set of nodes of this automaton is the product (computed using BDD product) state space */
            nodes_ = other->nodes_ * rabin->stateSpace_;
            /* the system nodes in the product are those nodes whose respective component in "other" belongs to the system nodes in "other" */
            sys_nodes_ = other->sys_nodes_;
            /* similarly, the environment nodes */
            env_nodes_ = other->env_nodes_;
            /* the set of variables is the union of variables of the arena and the rabin automaton */
            preVars_ = other->preVars_;
            for (auto i = rabin->stateVars_.begin(); i != rabin->stateVars_.end(); ++i) {
                preVars_.push_back(*i);
            }
            postVars_ = other->postVars_;
            for (auto i = rabin->postStateVars_.begin(); i != rabin->postStateVars_.end(); ++i) {
                postVars_.push_back(*i);
            }
            /* the joint transitions are the product of the two respective transitions */
            tr_ = other->tr_ * rabin->transitions_;
            live_ = other->live_ * rabin->transitions_;
            /* compute the helper BDDs for synthesis */
            cubePost_ = sylvan::Bdd::VariablesCube(sort(postVars_));
            cubeOther_ = sylvan::Bdd::VariablesCube(sort(other_variables));
            tr_domain_ = tr_.ExistAbstract(cubePost_ * cubeOther_);
            live_domain_ = live_.ExistAbstract(cubePost_ * cubeOther_);
            /* the rabin pairs are directly inherited from the rabin automaton */
            RabinPairs_ = rabin->RabinPairs_;
        }
        /* function: Rabin
    *
    *  computes the fair adversarial rabin winning domain
    *
    * Input:
    *  accl_on   - true/false setting the accelerated fixpoint on/off
    *  M         - the bound on the iteration count for memorizing the BDDs from the past iterations
    *  initial_seed  - initial seed for warm starting the nu fixedpoints (for example the under-approximation fixpoint can be warm-started from the result of the over-approximation fixpoint)
    *  verbose   - the verbosity level (0-2, default=0)
    */
        sylvan::BDD Rabin(const bool accl_on,
                          const size_t M, /* the bound on the iteration count for memorizing the BDDs from the past iterations */
                          const sylvan::Bdd initial_seed,
                          const int verbose = 0);
        /* function: cpre
    *
    * computes the controllable predecessor
    *
    * cpre(Zi) = { (x,x') | (x,x') in tr_
    *                           &
    *                   x in sys_nodes_ => x' in Zi
    *                           &
    *                   x in env_nodes_ => forall x''. (x,x'') in tr_ => x'' in Zi
    *
    */
        sylvan::Bdd cpre(const sylvan::Bdd &Zi) {
            sylvan::Bdd Z = Zi;
            /* project onto state alphabet */
            Z = Z.ExistAbstract(cubePost_ * cubeOther_);
            /* swap variables */
            Z = Z.Permute(preVars_, postVars_);
            /* the controllable system edges */
            sylvan::Bdd W0 = tr_ & Z & sys_nodes_;
            /* the controllable environment edges */
            sylvan::Bdd nZ = !Z & nodes_;
            /* the environment nodes having an outgoing edge outside Zi */
            sylvan::Bdd F = tr_.AndAbstract(nZ, cubePost_ * cubeOther_) & env_nodes_;
            /* the other environment nodes are controllable */
            sylvan::Bdd nF = !F;
            sylvan::Bdd W1 = tr_ & nF & env_nodes_;
            /* return all the controllable edges */
            return (W0 | W1);
        }

        /* function: apre
    *
    * computes the almost sure predecessor
    *
    * apre(Y,Z) = cpre(Y) & { (x,x') | x in live_domain_ & x' in Z}
    *                               OR
    *               cpre(Z)
    *
    */

        sylvan::Bdd apre(const sylvan::Bdd Y,
                         const sylvan::Bdd Z) {
            sylvan::Bdd Z2 = Z;
            /* project onto state alphabet */
            Z2 = Z2.ExistAbstract(cubePost_ * cubeOther_);
            /* swap variables */
            Z2 = Z2.Permute(preVars_, postVars_);
            /* the nodes from which there are live edges to Z */
            sylvan::Bdd W0 = live_.AndAbstract(Z2, cubePost_ * cubeOther_);
            /* remove from W0 those env vertices which have some transition outside Y */
            sylvan::Bdd P = cpre(Y);
            sylvan::Bdd W1 = W0 & P;
            /* the edges in cpre(Z) are also in apre(Y,Z) */
            sylvan::Bdd W2 = W1 | cpre(Z);
            return W2;
        }
        /* function: printTabs
    *  used to print n number of tabs on the standard I/O during printing the results
    */
        void printTabs(const int n) {
            for (int i = 0; i < n; i++) {
                std::cout << "\t";
            }
        }
        /* function: to_dec
    *  convert a number with base "base" to a decimal number
    * position 0 is MSB
    */

        size_t to_dec(const size_t base, const std::vector<size_t> number) {
            size_t N = 0;
            for (size_t i = 1; i <= number.size(); i++) {
                N += number[i - 1] * pow(base, number.size() - i);
            }
            return N;
        }

        /* function: factorial
    *  compute the factorial
    */
        size_t factorial(const size_t n) {
            size_t ans = 1;
            for (int i = n; i > 1; i--)
                ans *= i;
            return ans;
        }

        /* function: pad_zeros
    *  pad zeros to a vector "vec1" to make its size equal to n
    */
        inline std::vector<size_t> pad_zeros(const std::vector<size_t> vec1, const size_t n) {
            std::vector<size_t> vec2 = vec1;
            for (size_t i = 0; i < n; i++) {
                vec2.push_back(0);
            }
            return vec2;
        }
        /* function: check_threshold
    *  returns true if all the elements in the vector called "vec" are below a given threshold
    */
        inline bool check_threshold(const std::vector<size_t> &vec, const size_t th) {
            for (size_t i = 0; i < vec.size(); i++) {
                if (vec[i] > th) {
                    return false;
                }
            }
            return true;
        }

        inline size_t first_index(const std::vector<size_t> &vec, const size_t th) {
            size_t i = 0;
            while (vec[i] <= th)
                i++;
            return i;
        }
    }; /* close class def */

    TASK_DECL_4(sylvan::BDD, RabinRecurse, fairsyn::FixedPoint *, sylvan::Bdd *, struct const_arg_recursive_rabin *, struct nconst_arg_recursive_rabin *)
#define RabinRecurse(fp, controller, arg_const, arg_nconst) CALL(RabinRecurse, (fp), (controller), (arg_const), (arg_nconst))

    sylvan::BDD FixedPoint::Rabin(const bool accl_on,
                                  const size_t M, /* the bound on the iteration count for memorizing the BDDs from the past iterations */
                                  const sylvan::Bdd initial_seed,
                                  const int verbose) {
        /* copy the rabin pairs */
        std::vector<rabin_pair_> pairs = RabinPairs_;
        size_t nrp = pairs.size(); /* number of rabin pairs */
        /* initialize a pair of trivial bdd-s */
        sylvan::Bdd top = initial_seed;
        sylvan::Bdd bot = sylvan::Bdd::bddZero();
        /* the exact scheme as per the piterman paper */
        std::vector<std::map<std::vector<size_t>, sylvan::Bdd>> *hist_Y = new std::vector<std::map<std::vector<size_t>, sylvan::Bdd>>;
        std::vector<std::map<std::vector<size_t>, sylvan::Bdd>> *hist_X = new std::vector<std::map<std::vector<size_t>, sylvan::Bdd>>;
        std::map<std::vector<size_t>, size_t> *permutation = new std::map<std::vector<size_t>, size_t>;
        // /* debug */
        // TicToc timer;
        // timer.tic();
        // /* debug ends */
        if (accl_on) {
            /* if acceleration is on, then populate hist_Y and hist_X with the respective initial values */
            /* the FIRST index is related to the depth, which is at most nrp+1 (the outer layer does not participate in the caching operation) */

            *permutation = map_all_permutations(nrp);

            for (int i = 0; i < (*permutation).size(); i++) {

                std::map<std::vector<size_t>, sylvan::Bdd> temp;
                hist_Y->push_back(temp);
                hist_X->push_back(temp);
            }
        }

        // /* debug */
        // double time = timer.toc();
        // std::fstream logfile;
        // logfile.open("/local/POPL_experiments/varying_M/brain04.txt", std::fstream::app);
        // logfile << "Init_time=" << time << " ";
        // /* debug ends */

        /* create variables for remembering the current indices of the fixpoint variables and the indices of the rabin pairs */
        std::vector<size_t> *indexY = new std::vector<size_t>;
        std::vector<size_t> *indexX = new std::vector<size_t>;
        std::vector<size_t> *indexRP = new std::vector<size_t>;
        /* the controller */
        sylvan::Bdd C = sylvan::Bdd::bddZero();
        /* initialize the sets for the nu fixed point */
        sylvan::Bdd Y = sylvan::Bdd::bddZero();
        sylvan::Bdd YY = initial_seed;
        for (int i = 0; Y.ExistAbstract(cubePost_ * cubeOther_) != YY.ExistAbstract(cubePost_ * cubeOther_); i++) {
            Y = YY;

            if (accl_on)
                indexY->push_back(i);

            if (verbose == 2) {
                cout << std::endl;
                std::cout << "\t Y0, iteration " << i << ", states = ";
                print_bdd_info(Y, preVars_);
            } else if (verbose == 1) {
                sylvan::BddSet preBddVars = sylvan::BddSet::fromVector(preVars_);
                std::cout << "Y0, iteration " << i << ", states = " << Y.ExistAbstract(cubePost_ * cubeOther_).SatCount(preBddVars) << "\n";
            }
            /* reset the controller */
            C = sylvan::Bdd::bddZero();
            /* initialize the sets for the mu fixed point */
            sylvan::Bdd X = sylvan::Bdd::bddOne();
            sylvan::Bdd XX = sylvan::Bdd::bddZero();
            for (int k = 0; X.ExistAbstract(cubePost_ * cubeOther_) != XX.ExistAbstract(cubePost_ * cubeOther_); k++) {
                X = XX;

                if (accl_on)
                    indexX->push_back(k);

                if (verbose == 2) {
                    cout << std::endl;
                    std::cout << "\t X0, iteration " << k << ", states = ";
                    print_bdd_info(X, preVars_);
                } else if (verbose == 1) {
                    sylvan::BddSet preBddVars = sylvan::BddSet::fromVector(preVars_);
                    std::cout << "\t X0, iteration " << k << ", states = " << X.ExistAbstract(cubePost_ * cubeOther_).SatCount(preBddVars) << "\n";
                }

                sylvan::Bdd term = apre(Y, X);
                /* the state-input pairs added by the outermost loop get the smallest rank */
                sylvan::Bdd N = term & (!(C.ExistAbstract(cubePost_ * cubeOther_)));
                C |= N;
                /* recursively solve the rest of the fp */
                const_arg_recursive_rabin arg_const_new = {
                        accl_on,
                        M, /* the bound on the iteration count for memorizing the BDDs from the past iterations */
                        1, /* depth of the next recursion level */
                        pairs,
                        initial_seed,
                        verbose};
                nconst_arg_recursive_rabin arg_nconst_new = {sylvan::Bdd::bddOne(),
                                                             term.ExistAbstract(cubePost_ * cubeOther_),
                                                             indexRP,
                                                             indexY,
                                                             indexX,
                                                             permutation,
                                                             hist_Y,
                                                             hist_X};
                XX = RUN(RabinRecurse, this, &C, &arg_const_new, &arg_nconst_new);
                if (accl_on)
                    indexX->pop_back();
            }
            YY = X;
            if (accl_on)
                indexY->pop_back();
        }
        /* the winning strategy is the set of controlled edges whose source belong to the system */

        if (verbose == 2) {
            print_bdd_info(C, preVars_, postVars_);
            std::cout << std::endl;
            std::cout << "\t sys_nodes_ , states = ";
            uint32_t nvars = (preVars_.size());
            sylvan::BddSet preBddVars = sylvan::BddSet::fromVector(preVars_);
            for (int i = 0; i < pow(2, preVars_.size()); i++) {
                sylvan::Bdd cur_node = ElementToBdd(i, preBddVars, nvars);
                sylvan::Bdd intermediate = sys_nodes_ & cur_node;
                if (cur_node == intermediate) {
                    cout << i << " ";
                }
            }
            cout << -1 << " ";
            cout << std::endl;
        }

        sylvan::Bdd strategy = (C & sys_nodes_);

        return strategy.GetBDD();
    }

    TASK_DECL_5(sylvan::BDD, RabinRecurseForLoop, size_t, fairsyn::FixedPoint *, sylvan::Bdd *, struct const_arg_recursive_rabin *, struct nconst_arg_recursive_rabin *);
#define RabinRecurseForLoop(i, fp, controller, arg_const, arg_nconst) (CALL(RabinRecurseForLoop, (i), (fp), (controller), (arg_const), (arg_nconst)))

    TASK_IMPL_4(sylvan::BDD,
                RabinRecurse,
                fairsyn::FixedPoint *, fp,
                sylvan::Bdd *, controller,
                struct const_arg_recursive_rabin *, arg_const,
                struct nconst_arg_recursive_rabin *, arg_nconst) {
        /* initialize a vector for storing the winning domain computed in each thread */

        std::vector<sylvan::Bdd> Y;
        std::vector<sylvan::Bdd *> C;
        for (size_t i = 0; i < arg_const->pairs.size(); i++) {
            Y.push_back(sylvan::Bdd::bddZero());
            sylvan::Bdd controller_copy = *controller;
            C.push_back(&controller_copy);
        }

        /* spawn as many additional worker threads as the number of remaining pairs minus 1 */
        for (size_t i = 1; i < arg_const->pairs.size(); i++) sylvan::bdd_refs_spawn(SPAWN(RabinRecurseForLoop, i, fp, C[i], arg_const, arg_nconst));

        /* the current worker thread is repsonsible for the first rabin pair */
        sylvan::BDD y = CALL(RabinRecurseForLoop, 0, fp, C[0], arg_const, arg_nconst);
        Y[0] = y;
        sylvan::bdd_refs_push(y);

        /* synchronize all the additional threads */
        for (size_t i = arg_const->pairs.size() - 1; i >= 1; i--) Y[i] = sylvan::bdd_refs_sync(SYNC(RabinRecurseForLoop));

        /* dereference the refs pushed during the call operation */
        sylvan::bdd_refs_pop(1);

        /* compute the union of all the winning domains and the controllers computed by the different worker threads */
        sylvan::Bdd U = sylvan::Bdd::bddZero();
        for (size_t i = 0; i < arg_const->pairs.size(); i++) {
            U |= Y[i];
            sylvan::Bdd N = *C[i] & (!(controller->ExistAbstract(fp->cubePost_ * fp->cubeOther_)));
            *controller |= N;
        }

        /* return the underlying MTBDD of U */
        return U.GetBDD();
    }

    ///* Sequential rabin fixpoint */

    TASK_IMPL_5(sylvan::BDD,
                RabinRecurseForLoop,
                const size_t, i,
                fairsyn::FixedPoint *, fp,
                sylvan::Bdd *, controller,
                struct const_arg_recursive_rabin *, arg_const,
                struct nconst_arg_recursive_rabin *, arg_nconst) {
        /* unpack the inputs */
        const bool accl_on = arg_const->accl_on;
        const size_t M = arg_const->M; /* the bound on the iteration count for memorizing the BDDs from the past iterations */
        const int depth = arg_const->depth;
        const std::vector<rabin_pair_> pairs = arg_const->pairs;
        const sylvan::Bdd initial_seed = arg_const->initial_seed;
        sylvan::Bdd seqR = arg_nconst->seqR & fp->tr_;
        sylvan::Bdd right = arg_nconst->right & fp->tr_;
        /* the original scheme from piterman pnueli paper */
        std::map<std::vector<size_t>, size_t> *permutation = arg_nconst->permutation;
        std::vector<std::map<std::vector<size_t>, sylvan::Bdd>> *hist_Y = arg_nconst->hist_Y;
        std::vector<std::map<std::vector<size_t>, sylvan::Bdd>> *hist_X = arg_nconst->hist_X;
        /* the scheme used in the Mascot-SDS paper */

        /* create variables for remembering the current indices of the fixpoint variables and the indices of the rabin pairs */
        std::vector<size_t> *indexY = new std::vector<size_t>;
        *indexY = *(arg_nconst->indexY);
        std::vector<size_t> *indexX = new std::vector<size_t>;
        *indexX = *(arg_nconst->indexX);
        std::vector<size_t> *indexRP = new std::vector<size_t>;
        *indexRP = *(arg_nconst->indexRP);

        if (accl_on)
            indexRP->push_back(pairs[i].rabin_index_);

        sylvan::Bdd G = pairs[i].G_;
        sylvan::Bdd nR = pairs[i].nR_;
        std::vector<rabin_pair_> remPairs = pairs;
        remPairs.erase(remPairs.begin() + i);

        const int verbose = arg_const->verbose;
        if (verbose >= 2) {
            fp->printTabs(3 * depth - 1);
            std::cout << "\n Current rabin pair: {";
            print_bdd_info(G, fp->preVars_);
            std::cout << "}, {";
            print_bdd_info(!nR & fp->nodes_, fp->preVars_);
            std::cout << "}\n";
        }
        if (verbose >= 1) {
            fp->printTabs(3 * depth - 1);
            std::cout << "Remaining pairs " << pairs.size() << "\n\n";
        }
        /* initialize a local copy for the controller */
        sylvan::Bdd C = sylvan::Bdd::bddZero();

        /////////////////////////  INITIALISING BDDS in YY

        /* initialize the sets for the nu fixed point */
        sylvan::Bdd Y = sylvan::Bdd::bddZero();
        sylvan::Bdd YY;
        if (accl_on) {

            std::vector<size_t> temp;
            std::vector<size_t> temp_perm;
            int cur = 0;

            while ((*indexX)[cur] <= M - 1 && cur < (*indexX).size()) {
                temp.push_back((*indexX)[cur]);
                temp_perm.push_back((*indexRP)[cur]);
                cur++;
            }

            if (cur) {

                size_t rank = (*permutation)[temp_perm];

                if ((*hist_Y)[rank].find(temp) != (*hist_Y)[rank].end()) {

                    YY = (*hist_Y)[rank][temp];

                } else {

                    (*hist_Y)[rank][temp] = initial_seed;

                    YY = initial_seed;
                }

            } else {
                YY = initial_seed;
            }

        } else {
            YY = initial_seed;
        }

        for (int j = 0; Y.ExistAbstract(fp->cubePost_ * fp->cubeOther_) != YY.ExistAbstract(fp->cubePost_ * fp->cubeOther_); j++) {

            Y = YY;

            if (accl_on)
                indexY->push_back(j);

            if (verbose == 2) {
                cout << std::endl;
                fp->printTabs(3 * depth);
                std::cout << "Y" << depth << ", iteration " << j << ", states = ";
                print_bdd_info(Y, fp->preVars_);
            } else if (verbose >= 1) {
                fp->printTabs(3 * depth);
                sylvan::BddSet preBddVars = sylvan::BddSet::fromVector(fp->preVars_);
                std::cout << "Y" << depth << ", iteration " << j << ", states = " << Y.ExistAbstract(fp->cubePost_ * fp->cubeOther_).SatCount(preBddVars) << "\n";
            }

            sylvan::Bdd term1 = (right | (seqR & (nR & (G & fp->cpre(Y)))));


            /* reset the local copy of the controller to the most recently added state-input pairs */
            sylvan::Bdd N = term1 & (!(controller->ExistAbstract(fp->cubePost_ * fp->cubeOther_)));
            C = *controller | N;

            /////////////////////////  INITIALISING BDDS in XX

            /* initialize the sets for the mu fixed point */
            sylvan::Bdd X = sylvan::Bdd::bddOne();
            sylvan::Bdd XX;

            if (accl_on && fp->check_threshold(*indexY, M - 1)) {

                size_t rank = (*permutation)[*indexRP];

                if ((*hist_X)[rank].find(*indexY) != (*hist_X)[rank].end()) {

                    XX = (*hist_X)[rank][*indexY];

                } else {

                    (*hist_X)[rank][*indexY] = sylvan::Bdd::bddZero();

                    XX = sylvan::Bdd::bddZero();
                }

            } else {

                XX = sylvan::Bdd::bddZero();
            }

            for (int k = 0; X.ExistAbstract(fp->cubePost_ * fp->cubeOther_) != XX.ExistAbstract(fp->cubePost_ * fp->cubeOther_); k++) {
                X = XX;
                if (accl_on)
                    indexX->push_back(k);
                if (verbose == 2) {
                    cout << std::endl;
                    fp->printTabs(3 * depth + 1);
                    std::cout << "X" << depth << ", iteration " << k << ", states = ";
                    print_bdd_info(X, fp->preVars_);
                } else if (verbose == 1) {
                    sylvan::BddSet preBddVars = sylvan::BddSet::fromVector(fp->preVars_);
                    fp->printTabs(3 * depth + 1);
                    std::cout << "X" << depth << ", iteration " << k << ", states = " << X.ExistAbstract(fp->cubePost_ * fp->cubeOther_).SatCount(preBddVars) << "\n";
                }
                sylvan::Bdd term2 = (term1 | (seqR & (nR & fp->apre(Y, X))));
                /* add the recently added state-input pairs to the controller */
                N = term2 & (!(C.ExistAbstract(fp->cubePost_ * fp->cubeOther_)));
                C |= N;
                if (remPairs.size() == 0) {
                    XX = term2;
                } else {
                    const_arg_recursive_rabin arg_const_new = {
                            accl_on,
                            M, /* the bound on the iteration count for memorizing the BDDs from the past iterations */
                            depth + 1,
                            remPairs,
                            initial_seed,
                            verbose};
                    nconst_arg_recursive_rabin arg_nconst_new = {(seqR & nR).ExistAbstract(fp->cubePost_ * fp->cubeOther_),
                                                                 term2.ExistAbstract(fp->cubePost_ * fp->cubeOther_),
                                                                 indexRP,
                                                                 indexY,
                                                                 indexX,
                                                                 permutation,
                                                                 hist_Y,
                                                                 hist_X};
                    XX = CALL(RabinRecurse, fp, &C, &arg_const_new, &arg_nconst_new);
                }
                if (accl_on)
                    indexX->pop_back();
            }

            YY = XX;

            /////////////////////////  UPDATING BDDS in XX

            if (accl_on) {

                size_t rank = (*permutation)[*indexRP];

                if (fp->check_threshold(*indexY, M - 1)) {

                    //                   assert ((*hist_X)[rank][*indexY]  <= (XX.ExistAbstract(fp->cubePost_ * fp->cubeOther_) * fp->tr_)) ;
                    if ((*hist_X)[rank][*indexY] <= (XX.ExistAbstract(fp->cubePost_ * fp->cubeOther_) * fp->tr_))
                        (*hist_X)[rank][*indexY] = (XX.ExistAbstract(fp->cubePost_ * fp->cubeOther_) * fp->tr_);
                }

                indexY->pop_back();
            }
        }

        *controller = C;

        /////////////////////////  UPDATING BDDS in YY

        if (accl_on) {

            size_t rank = (*permutation)[*indexRP];
            if (fp->check_threshold(*indexX, M - 1)) {
                //                assert((*hist_Y)[rank][*indexX]  >= (YY.ExistAbstract(fp->cubePost_ * fp->cubeOther_) * fp->tr_) );
                if ((*hist_Y)[rank][*indexX] >= (YY.ExistAbstract(fp->cubePost_ * fp->cubeOther_) * fp->tr_))
                    (*hist_Y)[rank][*indexX] = (YY.ExistAbstract(fp->cubePost_ * fp->cubeOther_) * fp->tr_);
            }
        }


        if (accl_on)
            indexRP->pop_back();

        YY = YY.ExistAbstract(fp->cubePost_ * fp->cubeOther_);

        return YY.GetBDD();
    }

}// namespace fairsyn
