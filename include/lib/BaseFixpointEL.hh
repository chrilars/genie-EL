/** @file BaseFixpoint.hh
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


namespace genie {
    template<class UBDD>
    struct const_arg_recursive_rabin {
        const bool accl_on;
        const size_t M; /* the bound on the iteration count for memorizing the BDDs from the past iterations */
        const int depth;
        const std::vector<genie::rabin_pair_<UBDD>> pairs;
        const UBDD initial_seed;
        const int verbose;
    };

    template<class UBDD>
    struct nconst_arg_recursive_rabin {
        UBDD seqR;
        UBDD right;
        std::vector<size_t> *indexRP;
        std::vector<size_t> *indexY;
        std::vector<size_t> *indexX;

        //////////////////// for acceleration
        std::vector<std::vector<std::vector<std::vector<UBDD>>>> *hist_Y;
        std::vector<std::vector<std::vector<std::vector<UBDD>>>> *hist_X;
    };

    /**
     * @brief base class for the fixed point computation for the Rabin specification
     *
     * provides the fixed point computation for the Rabin specification
     * on finite transition systems with the edge fairness condition
     */
    template<class UBDD>
    class BaseFixpoint {
    public:
        UBDD base_;                                 /**< the bdd manager */
        const char *str_;                           /**< mode of calculation */
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

        /**
         * @brief cube of all the variables
         * */
        virtual inline UBDD CubeNotState() = 0;

        /**
         * @brief computes the controllable predecessor
         * @details
         * cpre(Zi) = { (x,u) | exists w, exists x': (x,u,w,x') in transitionRelation
         *                    and forall w, forall x': (x,u,w,x') in transitionRelation  => x' in Zi }
         */
        virtual UBDD cpre(const UBDD &Zi) = 0;

        /**
         * @brief computes the almost sure predecessor
         * @details
         * apre(Y,Z) = { (x,u) | forall w, forall x': (x,u,w,x') in maybeTransition_ => x' in Y
         *                    and forall w, exists x' in Z: (x,u,w,x') in sureTransition_ } OR cpre(Z)
         */
        virtual UBDD apre(const UBDD &Y, const UBDD &Z) = 0;

        virtual void print_rabin_info(const UBDD &store,
                                      const char *mode,
                                      int verbose,
                                      int iteration = 0,
                                      int depth = 0) = 0;

        /**
         * @brief computes the fair adversarial rabin winning domain
         * @param accl_on   - true/false setting the accelerated fixpoint on/off
         * @param M         - the bound on the iteration count for memorizing the BDDs from the past iterations
         * @param initial_seed  - initial seed for warm starting the nu fixpoints (for example the under-approximation fixpoint can be warm-started from the result of the over-approximation fixpoint)
         * @param verbose   - the verbosity level (0-2, default=0)
         */
        UBDD Rabin(const bool accl_on,
                   const size_t M, /* the bound on the iteration count for memorizing the BDDs from the past iterations */
                   const UBDD &initial_seed,
                   const int verbose,
                   std::function<UBDD(BaseFixpoint<UBDD> *, UBDD &, const_arg_recursive_rabin<UBDD>,
                                      nconst_arg_recursive_rabin<UBDD>)> RR = SequentialRabinRecurse) {
            /* copy the rabin pairs */
            std::vector<rabin_pair_<UBDD>> pairs = RabinPairs_;
            size_t nrp = pairs.size(); /* number of rabin pairs */
            /* initialize a pair of trivial bdd-s */
            UBDD top = initial_seed;
            UBDD bot = base_.zero();
            /* the exact scheme as per the piterman paper */
            auto *hist_Y = new std::vector<std::vector<std::vector<std::vector<UBDD>>>>;
            auto *hist_X = new std::vector<std::vector<std::vector<std::vector<UBDD>>>>;

            /* create variables for remembering the current indices of the fixpoint variables and the indices of the rabin pairs */
            auto *indexY = new std::vector<size_t>;
            auto *indexX = new std::vector<size_t>;
            auto *indexRP = new std::vector<size_t>;
            /* the controller */
            UBDD C = base_.zero();
            /* initialize the sets for the nu fixed point */
            UBDD Y = base_.zero();
            UBDD YY = initial_seed;
            for (int i = 0; Y.existAbstract(CubeNotState()) != YY.existAbstract(CubeNotState()); i++) {
                Y = YY;

                print_rabin_info(Y, "Y", verbose, i);

                /* reset the controller */
                C = base_.zero();
                /* initialize the sets for the mu fixed point */
                UBDD X = base_.one();
                UBDD XX = base_.zero();
                for (int k = 0; X.existAbstract(CubeNotState()) != XX.existAbstract(CubeNotState()); k++) {
                    X = XX;


                    print_rabin_info(X, "X", verbose, k);

                    UBDD term;
                    term = apre(Y, X);
                    /* the state-input pairs added by the outermost loop get the smallest rank */
                    UBDD N = term & (!(C.existAbstract(CubeNotState())));
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
                                                                       term.existAbstract(CubeNotState()),
                                                                       indexRP,
                                                                       indexY,
                                                                       indexX,
                                                                       hist_Y,
                                                                       hist_X};
                    //                    XX = RUN(RabinRecurse, this, &C, &arg_const_new, &arg_nconst_new);
                    XX = RR(this, C, arg_const_new, arg_nconst_new);
                }
                YY = X;
            }
            /* the winning strategy is the set of controlled edges whose source belong to the system */

            print_rabin_info(C, "end", verbose);
            return C;
        }


        /*
        ELSeq(term1,zielonka_tree,t (node in zielonka tree)){
            if t is winning:
                X = V
                XX = emptyset
            else:
                X = emptyset
                XX = V
            
            while X != XX:
                if leaf:
                    return leaf and term1
                else:
                    if t winning:
                        for s in R(t):
                            intersection recurse(term1,ZT,s)
                    else:
                        union recurse(term1,ZT,s)
        }
        */
        static UBDD SequentialRabinRecurse(BaseFixpoint<UBDD> *fp,
                                           UBDD &controller,
                                           const_arg_recursive_rabin<UBDD> rrConst,
                                           nconst_arg_recursive_rabin<UBDD> rrVars) {
            /* initialize the final solution to be returned in the end */
            /* unpack the inputs */
            const bool accl_on = rrConst.accl_on; // Unimportant
            const size_t M = rrConst.M; /* the bound on the iteration count for memorizing the BDDs from the past iterations */
            const int depth = rrConst.depth;
            const int verbose = rrConst.verbose;
            auto pairs = rrConst.pairs;
            auto initial_seed = rrConst.initial_seed;
            auto seqR = rrVars.seqR;
            auto right = rrVars.right;
            auto hist_Y = rrVars.hist_Y;
            auto hist_X = rrVars.hist_X;
            auto indexY = rrVars.indexY;
            auto indexX = rrVars.indexX;
            auto indexRP = rrVars.indexRP;

            // Temp variables as placeholders
            bool winning = false; // To be replaced by check on current Zielonka node

            UBDD U; // zero or one, depends on losing or winning, U == X_s
            if (winning):
                U = fp->base_.one();
            else:
                U = fp->base_.zero();

            // if R(s) != empty_set ? (And else will likely only contain return or some small stuff since cpre in this?)
            for (size_t i = 0; i < pairs.size(); i++) { // Assumes children because of the two simultaneous iterations, needs to be changed
                UBDD G = pairs[i].G_;
                UBDD nR = pairs[i].nR_;
                std::vector<rabin_pair_<UBDD>> remPairs = pairs;
                remPairs.erase(remPairs.begin() + i);
                /* initialize a local copy for the controller */
                UBDD C = fp->base_.zero();

                /* initialize the sets for the nu fixed point */ // if node in zielonka tree winning: else:
                UBDD Y = fp->base_.zero();
                UBDD YY = initial_seed;

                for (int j = 0; Y.existAbstract(fp->CubeNotState()) != YY.existAbstract(fp->CubeNotState()); j++) { // while X_s != W  (Keep only one of these)
                    Y = YY;
                    fp->print_rabin_info(Y, "Y", verbose, j, depth);

                    UBDD term1;
                    term1 = right | (seqR & nR & G & fp->cpre(Y)); // Do term1 calc later
                    /* reset the local copy of the controller to the most recently added state-input pairs */
                    UBDD N = term1 & (!(controller.existAbstract(fp->CubeNotState())));
                    C = controller | N;
                    /* initialize the sets for the mu fixed point */
                    UBDD X = fp->base_.one(); // Replace with U (== X_s) 
                    UBDD XX = fp->base_.zero(); // == 

                    for (int k = 0; X.existAbstract(fp->CubeNotState()) != XX.existAbstract(fp->CubeNotState()); k++) { // while X_s != W (Keep only one of these)
                        X = XX;

                        fp->print_rabin_info(X, "X", verbose, k, depth);
                        UBDD term2;
                        term2 = term1 | (seqR & nR & fp->apre(Y, X)); // Do we only need term1?
                        /* add the recently added state-input pairs to the controller */
                        N = term2 & (!(C.existAbstract(fp->CubeNotState()))); // For us (since no apre): term1 & ...
                        C |= N;
                        if (remPairs.empty()) { // if leaf in zielonka (if R(s) = empty_set)
                            XX = term2;
                        } else {
                            genie::const_arg_recursive_rabin<UBDD> arg_const_new = { // Most things in here not necessary
                                    accl_on,
                                    M, /* the bound on the iteration count for memorizing the BDDs from the past iterations */
                                    depth + 1,
                                    remPairs,
                                    initial_seed,
                                    verbose};
                            genie::nconst_arg_recursive_rabin<UBDD> arg_nconst_new = {
                                    seqR & nR,
                                    term2,
                                    indexRP,
                                    indexY,
                                    indexX,
                                    hist_Y,
                                    hist_X};
                            XX = SequentialRabinRecurse(fp, C, arg_const_new, arg_nconst_new); // need to pass along current node in zielonka
                        }
                    }
                    YY = XX;
                }
                if (winning): // Union or intersection depending on losing or winning
                    U &= YY;
                else:
                    U |= YY; 
                controller = C;
            }
            return U;
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
         *  @brief convert a number with base "base" to a decimal number position 0 is MSB
         */
        virtual size_t to_dec(const size_t base, const std::vector<size_t>& number) {
            size_t N = 0;
            for (size_t i = 1; i <= number.size(); i++) {
                //            N += number[i]*pow(base,i);
                N += number[i - 1] * pow(base, number.size() - i);
            }
            return N;
        }

        /**
         *  @brief compute the factorial
         */
        size_t factorial(const size_t n) {
            if (n == 0) {
                return 1;
            } else {
                return (n * factorial(n - 1));
            }
        }

        /**
         * @brief A utility function to count smaller characters on right of arr[low]
         */
        size_t findSmallerInRight(std::vector<size_t> str, size_t low) {
            size_t countRight = str[low]; /* number of characters smaller than the chosen one */
            for (size_t i = 0; i < low; i++) {
                if (str[i] < str[low]) {
                    countRight--;
                }
            }
            return countRight;
        }

        /**
         * @brief A function to find rank of a string in all permutations of characters
         */
        size_t findRank(size_t nsymbols, std::vector<size_t> str) {
            size_t len = str.size();
            size_t mul = factorial(len);
            size_t rank = 0;
            size_t countRight;
            for (size_t i = 0; i < len; ++i) {
                mul /= len - i;
                // count number of chars smaller than str[i] from str[i+1] to str[len-1]
                countRight = findSmallerInRight(str, i);
                rank += countRight * mul;
            }
            return rank;
        }

        /**
         * @brief check if all values are below threshold
         * @param vec - vector of elements
         * @param th  - threshold
         * @return true if all the elements in the vector called "vec" are below a given threshold
         */
        static inline bool check_threshold(const std::vector<size_t> &vec, const size_t th) {
            for (unsigned long i : vec) {
                if (i > th) {
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
                                             std::vector<std::vector<size_t>> &permutations,
                                             std::vector<size_t> &temp) {

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
}// namespace genie