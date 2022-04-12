/** @file BaseFixedPoint.hh
 */

#pragma once

#include "BaseRabinAutomaton.hh"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>
#include <sylvan_obj.hpp>
#include <vector>


namespace fairsyn {
    template <class UBDD>
    struct const_arg_recursive_rabin {
        const bool accl_on;
        const size_t M; /* the bound on the iteration count for memorizing the BDDs from the past iterations */
        const int depth;
        const std::vector<fairsyn::rabin_pair_<UBDD>> pairs;
        const UBDD initial_seed;
        const int verbose;
    };

    template <class UBDD>
    struct nconst_arg_recursive_rabin {
        UBDD seqR;
        UBDD right;
        std::vector<size_t> *indexRP;
        std::vector<size_t> *indexY;
        std::vector<size_t> *indexX;

        //////////////////// for acceleration
        std::map<std::vector<size_t>, size_t> *permutation;
        std::vector<std::map<std::vector<size_t>, UBDD>> *hist_Y;
        std::vector<std::map<std::vector<size_t>, UBDD>> *hist_X;
    };

    /**
     * @brief base class for the fixed point computation for the Rabin specification
     *
     * provides the fixed point computation for the Rabin specification
     * on finite transition systems with the edge fairness condition
     */
    template <class UBDD>
    class BaseFixedPoint {
    public:
        UBDD base_;                                 /**< the bdd manager */
        UBDD cubePost_;                             /**< cubes with post variables; used in the existential abstraction  */
        UBDD cubeOther_;                            /**< cubes with other variables (outside the pre and the post variables) on which the transitions possibly depend */
        std::vector<rabin_pair_<UBDD>> RabinPairs_; /**< vector of the rabin pairs */
        UBDD nodes_;                                /**< stores the nodes in a BDD */
        UBDD sys_nodes_;                            /**< stores the system nodes in a BDD */
        UBDD env_nodes_;                            /**< stores the environment nodes in a BDD */
        std::vector<size_t> preVars_;               /**< stores the "pre" node variable indices */
        std::vector<size_t> postVars_;              /**< stores the "post" node variable indices  */
        UBDD tr_;                                   /**< Transition BDD */
        UBDD live_;                                 /**< the live transitions (subset of the transition relation) */

        virtual UBDD RabinRecurse(UBDD controller,
                                  const_arg_recursive_rabin<UBDD> rrConst,
                                  nconst_arg_recursive_rabin<UBDD> rrVars) = 0;

        virtual void print_bdd_info(const UBDD &store,
                                    const std::vector<size_t> &preVars_,
                                    const int verbose) = 0;

        virtual void print_bdd_info(const UBDD &store,
                                    const std::vector<size_t> &preVars_,
                                    const std::vector<size_t> &postVars_,
                                    const int verbose) = 0;
        /**
         * @brief computes the controllable predecessor
         * @details
         * cpre(Zi) = { (x,u) | exists w, exists x': (x,u,w,x') in transitionRelation
         *                    and forall w, forall x': (x,u,w,x') in transitionRelation  => x' in Zi }
         */
        UBDD cpre(const UBDD &Zi) {
            UBDD Z = Zi;
            /* project onto state alphabet */
            Z = Z.existAbstract(cubePost_ * cubeOther_);
            /* swap variables */
            Z = Z.permute(preVars_, postVars_);
            /* the controllable system edges */
            UBDD W0 = tr_ & Z & sys_nodes_;
            /* the controllable environment edges */
            UBDD nZ = !Z & nodes_;
            /* the environment nodes having an outgoing edge outside Zi */
            UBDD F = tr_.andAbstract(nZ, cubePost_ * cubeOther_) & env_nodes_;
            /* the other environment nodes are controllable */
            UBDD nF = !F;
            UBDD W1 = tr_ & nF & env_nodes_;
            /* return all the controllable edges */
            return (W0 | W1);
        }

        /**
         * @brief computes the almost sure predecessor
         * @details
         * apre(Y,Z) = { (x,u) | forall w, forall x': (x,u,w,x') in maybeTransition_ => x' in Y
         *                    and forall w, exists x' in Z: (x,u,w,x') in sureTransition_ } OR cpre("maybe",Z)
         */
        UBDD apre(const UBDD &Y,
                  const UBDD &Z) {
            UBDD Z2 = Z;
            /* project onto state alphabet */
            Z2 = Z2.existAbstract(cubePost_ * cubeOther_);
            /* swap variables */
            Z2 = Z2.permute(preVars_, postVars_);
            /* the nodes from which there are live edges to Z */
            UBDD W0 = live_.andAbstract(Z2, cubePost_ * cubeOther_);
            /* remove from W0 those env vertices which have some transition outside Y */
            UBDD P = cpre(Y);
            UBDD W1 = W0 & P;
            /* the edges in cpre(Z) are also in apre(Y,Z) */
            UBDD W2 = W1 | cpre(Z);
            return W2;
        }

        /**
         * @brief computes the fair adversarial rabin winning domain
         * @param accl_on   - true/false setting the accelerated fixpoint on/off
         * @param M         - the bound on the iteration count for memorizing the BDDs from the past iterations
         * @param initial_seed  - initial seed for warm starting the nu fixedpoints (for example the under-approximation fixpoint can be warm-started from the result of the over-approximation fixpoint)
         * @param verbose   - the verbosity level (0-2, default=0)
         */
        UBDD Rabin(const bool accl_on,
                   const size_t M, /* the bound on the iteration count for memorizing the BDDs from the past iterations */
                   const UBDD &initial_seed,
                   const int verbose) {
            /* copy the rabin pairs */
            std::vector<rabin_pair_<UBDD>> pairs = RabinPairs_;
            size_t nrp = pairs.size(); /* number of rabin pairs */
            /* initialize a pair of trivial bdd-s */
            UBDD top = initial_seed;
            UBDD bot = base_.zero();
            /* the exact scheme as per the piterman paper */
            std::vector<std::map<std::vector<size_t>, UBDD>> *hist_Y = new std::vector<std::map<std::vector<size_t>, UBDD>>;
            std::vector<std::map<std::vector<size_t>, UBDD>> *hist_X = new std::vector<std::map<std::vector<size_t>, UBDD>>;
            std::map<std::vector<size_t>, size_t> *permutation = new std::map<std::vector<size_t>, size_t>;

            if (accl_on) {
                /* if acceleration is on, then populate hist_Y and hist_X with the respective initial values */
                /* the FIRST index is related to the depth, which is at most nrp+1 (the outer layer does not participate in the caching operation) */

                *permutation = map_all_permutations(nrp);

                for (int i = 0; i < (*permutation).size(); i++) {

                    std::map<std::vector<size_t>, UBDD> temp;
                    hist_Y->push_back(temp);
                    hist_X->push_back(temp);
                }
            }

            /* create variables for remembering the current indices of the fixpoint variables and the indices of the rabin pairs */
            std::vector<size_t> *indexY = new std::vector<size_t>;
            std::vector<size_t> *indexX = new std::vector<size_t>;
            std::vector<size_t> *indexRP = new std::vector<size_t>;
            /* the controller */
            UBDD C = base_.zero();
            /* initialize the sets for the nu fixed point */
            UBDD Y = base_.zero();
            UBDD YY = initial_seed;
            for (int i = 0; Y.existAbstract(cubePost_ * cubeOther_) != YY.existAbstract(cubePost_ * cubeOther_); i++) {
                Y = YY;

                if (accl_on)
                    indexY->push_back(i);

                if (verbose >= 1) {
                    std::cout << std::endl;
                    std::cout << "\t Y0, iteration " << i << ", states = ";
                    print_bdd_info(Y, preVars_, verbose);
                }

                /* reset the controller */
                C = base_.zero();
                /* initialize the sets for the mu fixed point */
                UBDD X = base_.one();
                UBDD XX = base_.zero();
                for (int k = 0; X.existAbstract(cubePost_ * cubeOther_) != XX.existAbstract(cubePost_ * cubeOther_); k++) {
                    X = XX;

                    if (accl_on)
                        indexX->push_back(k);

                    if (verbose >= 1) {
                        std::cout << std::endl;
                        std::cout << "\t X0, iteration " << k << ", states = ";
                        print_bdd_info(X, preVars_, verbose);
                    }

                    UBDD term = apre(Y, X);
                    /* the state-input pairs added by the outermost loop get the smallest rank */
                    UBDD N = term & (!(C.existAbstract(cubePost_ * cubeOther_)));
                    C |= N;
                    /* recursively solve the rest of the fp */
                    const_arg_recursive_rabin<UBDD> arg_const_new = {
                            accl_on,
                            M, /* the bound on the iteration count for memorizing the BDDs from the past iterations */
                            1, /* depth of the next recursion level */
                            pairs,
                            initial_seed,
                            verbose};
                    nconst_arg_recursive_rabin<UBDD> arg_nconst_new = {base_.one(),
                                                                       term.existAbstract(cubePost_ * cubeOther_),
                                                                       indexRP,
                                                                       indexY,
                                                                       indexX,
                                                                       permutation,
                                                                       hist_Y,
                                                                       hist_X};
                    //                    XX = RUN(RabinRecurse, this, &C, &arg_const_new, &arg_nconst_new);
                    XX = RabinRecurse(C, arg_const_new, arg_nconst_new);
                    if (accl_on)
                        indexX->pop_back();
                }
                YY = X;
                if (accl_on)
                    indexY->pop_back();
            }
            /* the winning strategy is the set of controlled edges whose source belong to the system */

            if (verbose == 2) {
                print_bdd_info(C, preVars_, postVars_, verbose);
                // todo
                //                std::cout << std::endl;
                //                std::cout << "\t sys_nodes_ , states = ";
                //
                //                uint32_t nvars = (preVars_.size());
                //                sylvan::BddSet preBddVars = sylvan::BddSet::fromVector(to_uint32_t(preVars_));
                //                for (int i = 0; i < pow(2, preVars_.size()); i++) {
                //                    UBDD cur_node = elementToBdd(i, preBddVars, nvars);
                //                    UBDD intermediate = sys_nodes_ & cur_node;
                //                    if (cur_node == intermediate) {
                //                        cout << i << " ";
                //                    }
                //                }
                //                cout << -1 << " ";
                //                cout << std::endl;
            }

            UBDD strategy = (C & sys_nodes_);

            return strategy;
        }

        /**
         *  @brief used to print n number of tabs on the standard I/O during printing the results
         */
        static void printTabs(const int n) {
            for (int i = 0; i < n; i++) {
                std::cout << "\t";
            }
        }

        /**
         * @brief check if all values are below threshold
         * @param vec - vector of elements
         * @param th  - threshold
         * @return true if all the elements in the vector called "vec" are below a given threshold
         */
        static inline bool check_threshold(const std::vector<size_t> &vec, const size_t th) {
            for (size_t i = 0; i < vec.size(); i++) {
                if (vec[i] > th) {
                    return false;
                }
            }
            return true;
        }

        static std::map<std::vector<size_t>, size_t> map_all_permutations(int nrp) {

            std::map<std::vector<size_t>, size_t> m;

            std::vector<std::vector<size_t>> all_permutations;

            std::vector<size_t> v, temp;
            std::vector<bool> vis;
            for (int i = 0; i < nrp; i++) {
                v.push_back(i);
                vis.push_back(false);
            }

            produce_all_permutations(v, vis, all_permutations, temp);

            for (int i = 0; i < all_permutations.size(); i++)
                m[all_permutations[i]] = i;

            return m;
        }

        static void produce_all_permutations(std::vector<size_t> &v, std::vector<bool> &vis,
                                             std::vector<std::vector<size_t>> &permutations, std::vector<size_t> &temp) {

            permutations.push_back(temp);

            if (temp.size() == v.size())
                return;

            for (int i = 0; i < v.size(); i++) {

                if (vis[i])
                    continue;

                temp.push_back(v[i]);
                vis[i] = true;
                produce_all_permutations(v, vis, permutations, temp);
                temp.pop_back();
            }
        }
    }; /* close class def */
}// namespace fairsyn
