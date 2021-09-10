/*
 * FixedPoint.hh
 *
 */

#pragma once

#include <iostream>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <vector>
#include <array>
#include <sylvan_obj.hpp>

#include "RabinAutomaton.hh"
#include "lib/Arena.hh"
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
    std::vector<std::vector<std::vector<std::vector<sylvan::Bdd>>>>* hist_Y;
    std::vector<std::vector<std::vector<std::vector<sylvan::Bdd>>>>* hist_X;
//    std::vector<std::vector<std::vector<sylvan::Bdd>>>* hist_Y;
//    std::vector<std::vector<std::vector<sylvan::Bdd>>>* hist_X;
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
class FixedPoint: public Arena {
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
    FixedPoint(std::vector<size_t>& nodes,
               std::vector<size_t>& sys_nodes,
               std::vector<size_t>& env_nodes,
               std::vector<std::vector<size_t>>& transitions,
               std::vector<std::vector<size_t>>& live_edges,
               std::vector<std::array<std::vector<size_t>,2>>& RabinPairs , std::set<uint32_t>& all_variables )  : Arena( all_variables, nodes, sys_nodes, env_nodes, transitions, live_edges)
    {
        cubePost_=sylvan::Bdd::VariablesCube(sort(postVars_));
        // std::vector<uint32_t> other_variables;
        // for (size_t i = 0; i < all_variables.size(); i++) {
        //     bool outside_support=true;
        //     for (size_t j = 0; j < preVars_.size(); j++) {
        //         if (all_variables[i]==preVars_[j]) {
        //             outside_support=false;
        //             break;
        //         }
        //     }
        //     if (outside_support) {
        //         for (size_t j = 0; j < postVars_.size(); j++) {
        //             if (all_variables[i]==postVars_[j]) {
        //                 outside_support=false;
        //                 break;
        //             }
        //         }
        //     }
        //     if (outside_support)
        //         other_variables.push_back(i);
        // }
        cubeOther_=sylvan::Bdd::bddOne();
        tr_domain_=tr_.ExistAbstract(cubePost_);
        live_domain_=live_.ExistAbstract(cubePost_);
        sylvan::BddSet preBddVars=sylvan::BddSet::fromVector(preVars_);
        for (size_t i=0; i<RabinPairs.size(); ++i) {
            rabin_pair_ p;
            p.rabin_index_=i;
            p.G_=SetToBdd(RabinPairs[i][0], preBddVars, preVars_.size());
            p.nR_=(!(SetToBdd(RabinPairs[i][1], preBddVars, preVars_.size())))&nodes_;
            RabinPairs_.push_back(p);
        }
    }
    /*
     * constructor
     */
    FixedPoint(const Arena* other,
               const RabinAutomaton* rabin) {
        /* the set of nodes of this automaton is the product (computed using BDD product) state space */
        nodes_=other->nodes_*rabin->stateSpace_;
        /* the system nodes in the product are those nodes whose respective component in "other" belongs to the system nodes in "other" */
        sys_nodes_=other->sys_nodes_;
        /* similarly, the environment nodes */
        env_nodes_=other->env_nodes_;
        /* the set of variables is the union of variables of the arena and the rabin automaton */
        preVars_=other->preVars_;
        for (auto i=rabin->stateVars_.begin(); i!=rabin->stateVars_.end(); ++i) {
            preVars_.push_back(*i);
        }
        postVars_=other->postVars_;
        for (auto i=rabin->postStateVars_.begin(); i!=rabin->postStateVars_.end(); ++i) {
            postVars_.push_back(*i);
        }
        /* the joint transitions are the product of the two respective transitions */
        tr_=other->tr_*rabin->transitions_;
        live_=other->live_*rabin->transitions_;
        /* compute the helper BDDs for synthesis */
        cubePost_=sylvan::Bdd::VariablesCube(sort(postVars_));
        cubeOther_=sylvan::Bdd::bddOne();
        tr_domain_=tr_.ExistAbstract(cubePost_);
        live_domain_=live_.ExistAbstract(cubePost_);
        /* the rabin pairs are directly inherited from the rabin automaton */
        RabinPairs_=rabin->RabinPairs_;
    }
    /*
     * constructor
     */
    FixedPoint(const Arena* other,
               const RabinAutomaton* rabin,
                const std::vector<uint32_t> other_variables) {
        /* the set of nodes of this automaton is the product (computed using BDD product) state space */
        nodes_=other->nodes_*rabin->stateSpace_;
        /* the system nodes in the product are those nodes whose respective component in "other" belongs to the system nodes in "other" */
        sys_nodes_=other->sys_nodes_;
        /* similarly, the environment nodes */
        env_nodes_=other->env_nodes_;
        /* the set of variables is the union of variables of the arena and the rabin automaton */
        preVars_=other->preVars_;
        for (auto i=rabin->stateVars_.begin(); i!=rabin->stateVars_.end(); ++i) {
            preVars_.push_back(*i);
        }
        postVars_=other->postVars_;
        for (auto i=rabin->postStateVars_.begin(); i!=rabin->postStateVars_.end(); ++i) {
            postVars_.push_back(*i);
        }
        /* the joint transitions are the product of the two respective transitions */
        tr_=other->tr_*rabin->transitions_;
        live_=other->live_*rabin->transitions_;
        /* compute the helper BDDs for synthesis */
        cubePost_=sylvan::Bdd::VariablesCube(sort(postVars_));
        cubeOther_=sylvan::Bdd::VariablesCube(sort(other_variables));
        tr_domain_=tr_.ExistAbstract(cubePost_*cubeOther_);
        live_domain_=live_.ExistAbstract(cubePost_*cubeOther_);
        /* the rabin pairs are directly inherited from the rabin automaton */
        RabinPairs_=rabin->RabinPairs_;
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
                           const int verbose=0);
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
    sylvan::Bdd cpre(const sylvan::Bdd& Zi)  {
        sylvan::Bdd Z = Zi;
        /* debug */
        // std::cout << "states = " << Z.ExistAbstract(cubePost_*cubeOther_).SatCount(preVars_.size()) << "\n";
        /* debug ends */
        /* project onto state alphabet */
        Z=Z.ExistAbstract(cubePost_*cubeOther_);
        /* debug */
        // std::cout << "states = " << Z.ExistAbstract(cubePost_*cubeOther_).SatCount(preVars_.size()) << "\n";
        /* debug ends */
        /* swap variables */
        Z=Z.Permute(preVars_, postVars_);
        /* the controllable system edges */
        sylvan::Bdd W0 = tr_ & Z & sys_nodes_;
        /* debug */
        // std::cout << "states = " << (tr_).ExistAbstract(cubePost_*cubeOther_).SatCount(preVars_.size()) << "\n";
        // sylvan::Bdd b = (tr_).ExistAbstract(cubePost_*cubeOther_).Support();
        // sylvan::Bdd cubePre = sylvan::Bdd::VariablesCube(sort(preVars_));
        // if ((b & (!cubePre)) != sylvan::Bdd::bddZero()) {
        //     std::cout << "Problem.\n";
        // }
        // std::vector<uint32_t> v = cubePost_.toVector();
        // std::vector<uint32_t> v2 = cubeOther_.toVector();
        /* debug ends */
        /* the controllable environment edges */
        sylvan::Bdd nZ = !Z & nodes_;
        /* the environment nodes having an outgoing edge outside Zi */
        sylvan::Bdd F = tr_.AndAbstract(nZ,cubePost_*cubeOther_) & env_nodes_;
        /* the other environment nodes are controllable */
        sylvan::Bdd nF = !F;
        sylvan::Bdd W1 = tr_ & nF & env_nodes_;
        /* debug */
        // std::cout << "states = " << W0.ExistAbstract(cubePost_*cubeOther_).SatCount(preVars_.size()) << "\n";
        // std::cout << "states = " << W1.ExistAbstract(cubePost_*cubeOther_).SatCount(preVars_.size()) << "\n";
        /* debug ends */
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
                     const sylvan::Bdd Z)  {
        sylvan::Bdd Z2 = Z;
        /* project onto state alphabet */
        Z2=Z2.ExistAbstract(cubePost_*cubeOther_);
        /* swap variables */
        Z2=Z2.Permute(preVars_, postVars_);
        /* the nodes from which there are live edges to Z */
        sylvan::Bdd W0 = live_.AndAbstract(Z2, cubePost_*cubeOther_);
        /* remove from W0 those env vertices which have some transition outside Y */
        sylvan::Bdd P = cpre(Y);
        sylvan::Bdd W1 = W0 & P;
        /* the edges in cpre(Z) are also in apre(Y,Z) */
        sylvan::Bdd W2 = W1 | cpre(Z);
        return W2;
    }
    /* function: printTabs
     *  used to print n number of tabs on the standard I/O during printing the results */
    void printTabs(const int n) {
        for (int i=0; i<n; i++) {
            std::cout << "\t";
        }
    }
    /* function: to_dec
     *  convert a number with base "base" to a decimal number
     * position 0 is MSB */
    size_t to_dec(const size_t base, const std::vector<size_t> number) {
        size_t N=0;
        for (size_t i=1; i<=number.size(); i++) {
//            N += number[i]*pow(base,i);
            N += number[i-1]*pow(base,number.size()-i);
        }
        return N;
    }
    /* function: factorial
     *  compute the factorial */
    size_t factorial(const size_t n) {
        if (n==0) {
            return 1;
        } else {
            return (n*factorial(n-1));
        }
    }
    /* function: lexi_order
     *  determine the position of a given permutation "seq1" in the lexicographically ordered set of all permutations of N numbers */
    size_t lexi_order(const std::vector<size_t> seq1, const size_t N) {
        /* first, make the sequence complete by filling the extra digits with the lexicographically smallest suffix */
        std::vector<size_t> seq2 = seq1;
        if (seq2.size()<N+1) {
            for (size_t i=0; i<=N; i++) {
                bool number_present=false;
                for (size_t j=0; j<seq2.size(); j++) {
                    if (seq2[j]==i) {
                        number_present=true;
                        break;
                    }
                }
                if (!number_present)
                    seq2.push_back(i);
            }
        }
        /* second, compute the Lehmer code of the resulting padded sequence */
        std::vector<size_t> lehmer;
        for (size_t i=0; i<seq2.size(); i++) {
            size_t right_inversion=0;
            for (size_t j=i+1; j<seq2.size(); j++) {
                if (seq2[i]>seq2[j]) {
                    right_inversion++;
                }
            }
            lehmer.push_back(right_inversion);
        }
        /* third, compute the decimal value of the Lehmer code, interpreting it as a factorial base number */
        size_t pos=0;
        for (size_t i=0; i<lehmer.size(); i++) {
            pos+= factorial(i)*lehmer[i];
        }
        return pos;
    }
    // A utility function to count smaller characters on right
    // of arr[low]
    size_t findSmallerInRight(std::vector<size_t> str, size_t low)
    {
        size_t countRight = str[low]; /* number of characters smaller than the chosen one */

        for (size_t i=0; i<low; i++) {
            if (str[i]<str[low]) {
                countRight--;
            }
        }

        return countRight;
    }

    // A function to find rank of a string in all permutations
    // of characters
    size_t findRank(size_t nsymbols, std::vector<size_t> str)
    {
        size_t len = str.size();
        size_t mul = factorial(len);
        size_t rank = 0;
        size_t countRight;

        for (size_t i = 0; i < len; ++i) {
            mul /= len - i;

            // count number of chars smaller than str[i]
            // from str[i+1] to str[len-1]
            countRight = findSmallerInRight(str, i);

            rank += countRight * mul;
        }

        return rank;
    }
    /* function: pad_zeros
     *  pad zeros to a vector "vec1" to make its size equal to n */
    inline std::vector<size_t> pad_zeros(const std::vector<size_t> vec1, const size_t n) {
        std::vector<size_t> vec2=vec1;
        for (size_t i=0; i<n; i++) {
            vec2.push_back(0);
        }
        return vec2;
    }
    /* function: check_threshold
     *  returns true if all the elements in the vector called "vec" are below a given threshold */
    inline bool check_threshold(const std::vector<size_t> vec, const size_t th) {
        for (size_t i=0; i<vec.size(); i++) {
            if (vec[i]>th) {
                return false;
            }
        }
        return true;
    }
}; /* close class def */

TASK_DECL_4(sylvan::BDD, RabinRecurse, fairsyn::FixedPoint*, sylvan::Bdd*, struct const_arg_recursive_rabin*, struct nconst_arg_recursive_rabin*)
#define RabinRecurse(fp, controller, arg_const, arg_nconst) CALL(RabinRecurse, (fp), (controller), (arg_const), (arg_nconst))

sylvan::BDD FixedPoint::Rabin(const bool accl_on,
                                const size_t M, /* the bound on the iteration count for memorizing the BDDs from the past iterations */
                                const sylvan::Bdd initial_seed,
                                const int verbose) {
    /* copy the rabin pairs */
    std::vector<rabin_pair_> pairs=RabinPairs_;
    size_t nrp = pairs.size(); /* number of rabin pairs */
    /* initialize a pair of trivial bdd-s */
    sylvan::Bdd top = initial_seed;
    sylvan::Bdd bot = sylvan::Bdd::bddZero();
    /* the exact scheme as per the piterman paper */
    std::vector<std::vector<std::vector<std::vector<sylvan::Bdd>>>>* hist_Y = new std::vector<std::vector<std::vector<std::vector<sylvan::Bdd>>>>;
    std::vector<std::vector<std::vector<std::vector<sylvan::Bdd>>>>* hist_X = new std::vector<std::vector<std::vector<std::vector<sylvan::Bdd>>>>;
    /* debug */
    TicToc timer;
    timer.tic();
    /* debug ends */
    if (accl_on) {
        /* if acceleration is on, then populate hist_Y and hist_X with the respective initial values */
        /* the FIRST index is related to the depth, which is at most nrp+1 (the outer layer does not participate in the caching operation) */
        for (size_t i=0; i<nrp; i++) {
            std::vector<std::vector<std::vector<sylvan::Bdd>>> x, y;
            size_t npos=i+1; /* the actual depth of the fixpoint variable, where the outermost variables have depth 0 */
            /* the SECOND index is the lexicographic position of the rabin pair subsequence 1...npos */
            for (size_t j=0; j<factorial(nrp)/factorial(nrp-npos); j++) {
//            for (size_t j=0; j<pow(nrp,npos); j++) {
                std::vector<std::vector<sylvan::Bdd>> xx, yy;
                /* the THIRD index is the value of the sequence 0...npos-1 of the respective variables (Y sequence for hist_Y and X sequence for hist_X) */
                for (size_t k=0; k<M; k++) {
                    std::vector<sylvan::Bdd> yyy(pow(M,npos), top); /* the FOURTH index for hist_Y is the value of the sequence 0...npos-1 of the X variables */
                    yy.push_back(yyy);
                    std::vector<sylvan::Bdd> xxx(pow(M,npos+1), bot); /* the FOURTH index for hist_X is the value of the sequence 0...npos of the Y variables */
                    xx.push_back(xxx);
                }
                y.push_back(yy);
                x.push_back(xx);
            }
            hist_Y->push_back(y);
            hist_X->push_back(x);
        }
    }
    /* debug */
    double time=timer.toc();
    std::fstream logfile;
    logfile.open("/local/POPL_experiments/varying_M/brain04.txt", std::fstream::app);
    logfile << "Init_time=" << time << " ";
    /* debug ends */
//    /* the scheme used in Mascot-SDS */
//    std::vector<std::vector<std::vector<sylvan::Bdd>>>* hist_Y = new std::vector<std::vector<std::vector<sylvan::Bdd>>>;
//    std::vector<std::vector<std::vector<sylvan::Bdd>>>* hist_X = new std::vector<std::vector<std::vector<sylvan::Bdd>>>;
//    if (accl_on) {
//        /* if acceleration is on, then populate hist_Y and hist_X with the respective initial values */
//        for (size_t i=0; i<factorial(nrp); i++) {
//            std::vector<std::vector<sylvan::Bdd>> x, y;
//            for (size_t j=0; j<pow(M,nrp+1); j++) {
//                std::vector<sylvan::Bdd> yy(nrp+1, top);
//                y.push_back(yy);
//                std::vector<sylvan::Bdd> xx(nrp+1, bot);
//                x.push_back(xx);
//            }
//            hist_Y->push_back(y);
//            hist_X->push_back(x);
//        }
//    }
//    /* New scheme in between Mascot-SDS and the original one */
//    std::vector<std::vector<std::vector<std::vector<sylvan::Bdd>>>>* hist_Y = new std::vector<std::vector<std::vector<std::vector<sylvan::Bdd>>>>;
//    std::vector<std::vector<std::vector<std::vector<sylvan::Bdd>>>>* hist_X = new std::vector<std::vector<std::vector<std::vector<sylvan::Bdd>>>>;
//    if (accl_on) {
//        /* if acceleration is on, then populate hist_Y and hist_X with the respective initial values */
//        for (size_t i=0; i<factorial(nrp); i++) {
//            std::vector<std::vector<std::vector<sylvan::Bdd>>> x, y;
//            for (size_t j=0; j<pow(M,nrp+1); j++) {
//                std::vector<std::vector<sylvan::Bdd>> xx, yy;
//                for (size_t k=0; k<nrp; k++) {
//                    std::vector<sylvan::Bdd> yyy(M, top);
//                    yy.push_back(yyy);
//                    std::vector<sylvan::Bdd> xxx(M, bot);
//                    xx.push_back(xxx);
//                }
//                x.push_back(xx);
//                y.push_back(yy);
//            }
//            hist_Y->push_back(y);
//            hist_X->push_back(x);
//        }
//    }

    /* create variables for remembering the current indices of the fixpoint variables and the indices of the rabin pairs */
    std::vector<size_t> *indexY = new std::vector<size_t>;
    std::vector<size_t> *indexX = new std::vector<size_t>;
    std::vector<size_t> *indexRP = new std::vector<size_t>;
    /* the controller */
    sylvan::Bdd C=sylvan::Bdd::bddZero();
    /* initialize the sets for the nu fixed point */
    sylvan::Bdd Y = sylvan::Bdd::bddZero();
    sylvan::Bdd YY = initial_seed;
    for (int i=0; Y.ExistAbstract(cubePost_*cubeOther_)!=YY.ExistAbstract(cubePost_*cubeOther_); i++) {
        Y=YY;
        if (accl_on)
            indexY->push_back(i);

            if (verbose == 2){
            		cout<<std::endl;
               std::cout << "\t Y0, iteration " << i << ", states = ";
               print_bdd_info(Y,preVars_);
               /* uint32_t nvars=preVars_.size();
                sylvan::BddSet preBddVars=sylvan::BddSet::fromVector(preVars_);
                for(int i=0;i<pow(2,preVars_.size());i++){
                	sylvan::Bdd cur_node=ElementToBdd(i,preBddVars,nvars);
                	sylvan::Bdd intermediate = cur_node & Y;
			    if( cur_node == intermediate ){
				cout<<i<<" ";
			    }
                }
                cout<<-1<<" ";
                cout<<std::endl;
            		*/
            }else if (verbose == 1) {
            // std::cout << "Y0, iteration " << i << ", states = " << Y.ExistAbstract(cubePost_).SatCount(preVars_.size()) << "\n";
            sylvan::BddSet preBddVars=sylvan::BddSet::fromVector(preVars_);
            std::cout << "Y0, iteration " << i << ", states = " << Y.ExistAbstract(cubePost_*cubeOther_).SatCount(preBddVars) << "\n";
        }
        /* reset the controller */
        C = sylvan::Bdd::bddZero();
        /* initialize the sets for the mu fixed point */
        sylvan::Bdd X = sylvan::Bdd::bddOne();
        sylvan::Bdd XX = sylvan::Bdd::bddZero();
        for (int k=0; X.ExistAbstract(cubePost_*cubeOther_)!=XX.ExistAbstract(cubePost_*cubeOther_); k++) {
            X=XX;
            if (accl_on)
                indexX->push_back(k);
            if (verbose == 2){
            		cout<<std::endl;
                std::cout << "\t X0, iteration " << k << ", states = " ;
                print_bdd_info(X,preVars_);
                /*
                uint32_t nvars=(preVars_.size());
                sylvan::BddSet preBddVars=sylvan::BddSet::fromVector(preVars_);
                for(int i=0;i<pow(2,preVars_.size());i++){
                	sylvan::Bdd cur_node=ElementToBdd(i,preBddVars,nvars);
                	sylvan::Bdd intermediate = cur_node & X;
			    if( cur_node == intermediate ){
				cout<<i<<" ";
			    }
                }
                cout<<-1<<" ";
                cout<<std::endl;*/

            }else if (verbose == 1) {
                // std::cout << "\t X0, iteration " << k << ", states = " << X.ExistAbstract(cubePost_).SatCount(preVars_.size()) << "\n";
                sylvan::BddSet preBddVars=sylvan::BddSet::fromVector(preVars_);
                std::cout << "\t X0, iteration " << k << ", states = " << X.ExistAbstract(cubePost_*cubeOther_).SatCount(preBddVars) << "\n";
            }
            sylvan::Bdd term = apre(Y,X);
            /* the state-input pairs added by the outermost loop get the smallest rank */
            sylvan::Bdd N = term & (!(C.ExistAbstract(cubePost_*cubeOther_)));
            C |= N;
            /* recursively solve the rest of the fp */
            const_arg_recursive_rabin arg_const_new = {
                accl_on,
                M, /* the bound on the iteration count for memorizing the BDDs from the past iterations */
                1, /* depth of the next recursion level */
                pairs,
                initial_seed,
                verbose
            };
            nconst_arg_recursive_rabin arg_nconst_new = {sylvan::Bdd::bddOne(),
                term.ExistAbstract(cubePost_*cubeOther_),
                indexRP,
                indexY,
                indexX,
                hist_Y,
                hist_X
            };
            XX= RUN(RabinRecurse, this, &C, &arg_const_new, &arg_nconst_new);
            if (accl_on)
                indexX->pop_back();
        }
        YY= X;
        if (accl_on)
            indexY->pop_back();
    }
    /* the winning strategy is the set of controlled edges whose source belong to the system */


    if (verbose == 2){
        print_bdd_info(C,preVars_,postVars_);
        std::cout<<std::endl;
        std::cout << "\t sys_nodes_ , states = ";
        uint32_t nvars=(preVars_.size());
        sylvan::BddSet preBddVars=sylvan::BddSet::fromVector(preVars_);
        for(int i=0;i<pow(2,preVars_.size());i++){
            sylvan::Bdd cur_node=ElementToBdd(i,preBddVars,nvars);
            sylvan::Bdd intermediate = sys_nodes_ & cur_node;
            if( cur_node == intermediate ) {
                cout<<i<<" ";
            }
        }
        cout<<-1<<" ";
        cout<<std::endl;
    }

    sylvan::Bdd strategy = (C & sys_nodes_);

    return strategy.GetBDD();
}

TASK_DECL_5(sylvan::BDD, RabinRecurseForLoop, size_t, fairsyn::FixedPoint*, sylvan::Bdd*, struct const_arg_recursive_rabin*, struct nconst_arg_recursive_rabin*);
#define RabinRecurseForLoop(i, fp, controller, arg_const, arg_nconst) (CALL(RabinRecurseForLoop, (i), (fp), (controller), (arg_const), (arg_nconst)))

TASK_IMPL_4(sylvan::BDD,
            RabinRecurse,
            fairsyn::FixedPoint*, fp,
            sylvan::Bdd*, controller,
            struct const_arg_recursive_rabin*, arg_const,
            struct nconst_arg_recursive_rabin*, arg_nconst) {
    /* initialize a vector for storing the winning domain computed in each thread */
    std::vector<sylvan::Bdd> Y;
    std::vector<sylvan::Bdd*> C;
    for (size_t i=0; i<arg_const->pairs.size(); i++) {
        Y.push_back(sylvan::Bdd::bddZero());
        sylvan::Bdd controller_copy=*controller;
        C.push_back(&controller_copy);
    }
    /* spawn as many additional worker threads as the number of remaining pairs minus 1 */
    for (size_t i=1; i<arg_const->pairs.size(); i++) sylvan::bdd_refs_spawn(SPAWN(RabinRecurseForLoop, i, fp, C[i], arg_const, arg_nconst));
    /* the current worker thread is repsonsible for the first rabin pair */
    sylvan::BDD y=CALL(RabinRecurseForLoop, 0, fp, C[0], arg_const, arg_nconst);
    Y[0]=y;
    sylvan::bdd_refs_push(y);
    /* synchronize all the additional threads */
    for (size_t i=arg_const->pairs.size()-1; i>=1; i--) Y[i]= sylvan::bdd_refs_sync(SYNC(RabinRecurseForLoop));
    /* dereference the refs pushed during the call operation */
    sylvan::bdd_refs_pop(1);
    /* compute the union of all the winning domains and the controllers computed by the different worker threads */
    sylvan::Bdd U=sylvan::Bdd::bddZero();
    for (size_t i=0; i<arg_const->pairs.size(); i++) {
        U|=Y[i];
        sylvan::Bdd N= *C[i] & (!(controller->ExistAbstract(fp->cubePost_*fp->cubeOther_)));
        *controller |= N;
    }
    /* return the underlying MTBDD of U */
    return U.GetBDD();
}

///* Sequential rabin fixpoint */
//TASK_IMPL_4(sylvan::BDD,
//            RabinRecurse,
//            fairsyn::FixedPoint*, fp,
//            sylvan::Bdd*, controller,
//            struct const_arg_recursive_rabin*, arg_const,
//            struct nconst_arg_recursive_rabin*, arg_nconst) {
//    /* initialize a vector for storing the winning domain computed in each thread */
//    std::vector<sylvan::Bdd> Y;
//    std::vector<sylvan::Bdd*> C;
//    for (size_t i=0; i<arg_const->pairs.size(); i++) {
//        Y.push_back(sylvan::Bdd::bddZero());
//        sylvan::Bdd controller_copy=*controller;
//        C.push_back(&controller_copy);
//    }
//    /* iterate over all the rabin pairs and recursively call the for loop while incrementally building the controller */
//    for (size_t i=0; i<arg_const->pairs.size(); i++) {
//        Y[i] = sylvan::bdd_refs_push(CALL(RabinRecurseForLoop, i, fp, C[i], arg_const, arg_nconst));
//    }
////    Y[0] = sylvan::bdd_refs_push(CALL(RabinRecurseForLoop, 0, fp, C[0], arg_const, arg_nconst));
//    /* compute the union of all the winning domains and the controllers computed by the different worker threads */
//    sylvan::Bdd U=sylvan::Bdd::bddZero();
//    for (size_t i=0; i<arg_const->pairs.size(); i++) {
//        U|=Y[i];
//        sylvan::Bdd N= *C[i] & (!(controller->ExistAbstract(fp->cubePost_)));
//        *controller |= N;
//    }
//    sylvan::bdd_refs_pop(arg_const->pairs.size());
//    /* return the underlying MTBDD of U */
//    return U.GetBDD();
//}

TASK_IMPL_5(sylvan::BDD,
            RabinRecurseForLoop,
            const size_t, i,
            fairsyn::FixedPoint*, fp,
            sylvan::Bdd*, controller,
            struct const_arg_recursive_rabin*, arg_const,
            struct nconst_arg_recursive_rabin*, arg_nconst) {
    /* unpack the inputs */
    const bool accl_on=arg_const->accl_on;
    const size_t M=arg_const->M; /* the bound on the iteration count for memorizing the BDDs from the past iterations */
    const int depth=arg_const->depth;
    const std::vector<rabin_pair_> pairs=arg_const->pairs;
    const sylvan::Bdd initial_seed=arg_const->initial_seed;
    sylvan::Bdd seqR=arg_nconst->seqR & fp->tr_;
    sylvan::Bdd right=arg_nconst->right & fp->tr_;
    /* the original scheme from piterman pnueli paper */
    std::vector<std::vector<std::vector<std::vector<sylvan::Bdd>>>>* hist_Y = arg_nconst->hist_Y;
    std::vector<std::vector<std::vector<std::vector<sylvan::Bdd>>>>* hist_X = arg_nconst->hist_X;
    /* the scheme used in the Mascot-SDS paper */
//    std::vector<std::vector<std::vector<sylvan::Bdd>>>* hist_Y = arg_nconst->hist_Y;
//    std::vector<std::vector<std::vector<sylvan::Bdd>>>* hist_X = arg_nconst->hist_X;
    /* create variables for remembering the current indices of the fixpoint variables and the indices of the rabin pairs */
    std::vector<size_t> *indexY = new std::vector<size_t>;
    *indexY=*(arg_nconst->indexY);
    std::vector<size_t> *indexX = new std::vector<size_t>;
    *indexX=*(arg_nconst->indexX);
    std::vector<size_t> *indexRP = new std::vector<size_t>;
    *indexRP=*(arg_nconst->indexRP);

    if (accl_on)
        indexRP->push_back(pairs[i].rabin_index_);
    sylvan::Bdd G = pairs[i].G_;
    sylvan::Bdd nR = pairs[i].nR_;
    std::vector<rabin_pair_> remPairs=pairs;
    remPairs.erase(remPairs.begin()+i);

    const int verbose=arg_const->verbose;
    if (verbose>=2) {
        fp->printTabs(3*depth-1);
        std::cout << "\n Current rabin pair: {";
        print_bdd_info(G, fp->preVars_);
        std::cout << "}, {";
        print_bdd_info(!nR & fp->nodes_, fp->preVars_);
        std::cout << "}\n";
    }
    if (verbose>=1) {
        fp->printTabs(3*depth-1);
        std::cout << "Remaining pairs " << pairs.size() << "\n\n";
    }
    /* initialize a local copy for the controller */
    sylvan::Bdd C = sylvan::Bdd::bddZero();;
    /* initialize the sets for the nu fixed point */
    sylvan::Bdd Y = sylvan::Bdd::bddZero();
    sylvan::Bdd YY;
    if (accl_on && fp->check_threshold(*indexY,M-1) && fp->check_threshold(*indexX,M-1)) {
//        YY = (*hist_Y)
//        [fp->lexi_order(*indexRP,fp->RabinPairs_.size()-1)]
//        [fp->to_dec(M,fp->pad_zeros(*indexX,remPairs.size()+1))]
//        [depth-1]
//        [std::min((*indexY)[0],M-1)];
//        YY = (*hist_Y)[fp->lexi_order(*indexRP,fp->RabinPairs_.size()-1)][fp->to_dec(M,fp->pad_zeros(*indexX,remPairs.size()+1))][depth-1];
        YY = (*hist_Y)[depth-1]\
        [fp->findRank(fp->RabinPairs_.size(),*indexRP)]\
        [std::min((*indexY)[0],M-1)]\
        [fp->to_dec(M,*indexX)];
    } else {
        YY = initial_seed;
    }
    for (int j=0; Y.ExistAbstract(fp->cubePost_*fp->cubeOther_)!=YY.ExistAbstract(fp->cubePost_*fp->cubeOther_); j++) {
        Y=YY;
        if (accl_on)
            indexY->push_back(j);
        if (verbose == 2){
            cout<<std::endl;
            fp->printTabs(3*depth);
            std::cout << "Y" << depth << ", iteration " << j << ", states = " ;
            print_bdd_info(Y,fp->preVars_);
        } else if (verbose>=1) {
            fp->printTabs(3*depth);
            // std::cout << "Y" << depth << ", iteration " << j << ", states = " << Y.ExistAbstract(fp->cubePost_*fp->cubeOther_).SatCount(fp->preVars_.size()) << "\n";
            sylvan::BddSet preBddVars=sylvan::BddSet::fromVector(fp->preVars_);
            std::cout << "Y" << depth << ", iteration " << j << ", states = " << Y.ExistAbstract(fp->cubePost_*fp->cubeOther_).SatCount(preBddVars) << "\n";
            // std::cout << "Y" << depth << ", iteration " << j << "\n";
        }
        sylvan::Bdd term1 = (right | (seqR & (nR & (G & fp->cpre(Y)))));
        // /* debug */
        // std::cout << "cpre = " << (fp->cpre(Y)).ExistAbstract(fp->cubePost_*fp->cubeOther_).SatCount(fp->preVars_.size()) << "\n";
        // std::cout << "states = " << (G & fp->cpre(Y)).ExistAbstract(fp->cubePost_*fp->cubeOther_).SatCount(fp->preVars_.size()) << "\n";
        // std::cout << "states = " << (nR & (G & fp->cpre(Y))).ExistAbstract(fp->cubePost_*fp->cubeOther_).SatCount(fp->preVars_.size()) << "\n";
        // std::cout << "states = " << (seqR & (nR & (G & fp->cpre(Y)))).ExistAbstract(fp->cubePost_*fp->cubeOther_).SatCount(fp->preVars_.size()) << "\n";
        // std::cout << "states = " << right.ExistAbstract(fp->cubePost_*fp->cubeOther_).SatCount(fp->preVars_.size()) << "\n";
        // std::cout << "states = " << term1.ExistAbstract(fp->cubePost_*fp->cubeOther_).SatCount(fp->preVars_.size()) << "\n";
        // /* debug ends */
        /* reset the local copy of the controller to the most recently added state-input pairs */
        sylvan::Bdd N = term1 & (!(controller->ExistAbstract(fp->cubePost_*fp->cubeOther_)));
        C=*controller | N;
        /* initialize the sets for the mu fixed point */
        sylvan::Bdd X = sylvan::Bdd::bddOne();
        sylvan::Bdd XX;
        // /* debug */
        // std::vector<size_t> v1={1,0,0};
        // std::vector<size_t> v2={0,2,1};
        // std::vector<size_t> v3={0,0,0,0};
        // if (*indexX==v1 && *indexRP==v2 && *indexY==v3) {
        //     cout << "Here";
        // }
        // /* debug end */
        if (accl_on && fp->check_threshold(*indexY,M-1) && fp->check_threshold(*indexX,M-1)) {
//            XX=(*hist_X)
//            [fp->lexi_order(*indexRP,fp->RabinPairs_.size()-1)]
//            [fp->to_dec(M,fp->pad_zeros(*indexY,remPairs.size()))]
//            [depth-1]
//            [std::min((*indexX)[0],M-1)];
//            XX=(*hist_X)[fp->lexi_order(*indexRP,fp->RabinPairs_.size()-1)][fp->to_dec(M,fp->pad_zeros(*indexY,remPairs.size()))][depth-1];
            XX=(*hist_X)[depth-1]\
            [fp->findRank(fp->RabinPairs_.size(),*indexRP)]\
            [std::min((*indexX)[0],M-1)]\
            [fp->to_dec(M,*indexY)];
        } else {
            XX = sylvan::Bdd::bddZero();
        }
        // for (int k=0; X!=XX; k++) {
        for (int k=0; X.ExistAbstract(fp->cubePost_*fp->cubeOther_)!=XX.ExistAbstract(fp->cubePost_*fp->cubeOther_); k++) {
            X=XX;
            if (accl_on)
                indexX->push_back(k);
            if (verbose == 2){
            		cout<<std::endl;
            	fp->printTabs(3*depth+1);
                std::cout << "X" << depth << ", iteration " << k << ", states = " ;
                print_bdd_info(X,fp->preVars_);
            } else if (verbose==1) {
                sylvan::BddSet preBddVars=sylvan::BddSet::fromVector(fp->preVars_);
                fp->printTabs(3*depth+1);
                // std::cout << "X" << depth << ", iteration " << k << ", states = " << X.ExistAbstract(fp->cubePost_*fp->cubeOther_).SatCount(fp->preVars_.size()) << "\n";
                // sylvan::Bdd test = X.ExistAbstract(fp->cubePost_*fp->cubeOther_);
                // std::cout << "X" << depth << ", iteration " << k << ", states = " << test.SatCount(preBddVars) << "\n";
                std::cout << "X" << depth << ", iteration " << k << ", states = " << X.ExistAbstract(fp->cubePost_*fp->cubeOther_).SatCount(preBddVars) << "\n";
                // std::cout << depth << k << "\n";
            }
            sylvan::Bdd term2 = (term1 | (seqR & (nR & fp->apre(Y,X))));
            /* debug */
            // std::cout << "apre = " << (fp->apre(Y,X)).ExistAbstract(fp->cubePost_*fp->cubeOther_).SatCount(fp->preVars_.size()) << "\n";
            // std::cout << "states = " << (nR & fp->apre(Y,X)).ExistAbstract(fp->cubePost_*fp->cubeOther_).SatCount(fp->preVars_.size()) << "\n";
            // std::cout << "states = " << (seqR & (nR & fp->apre(Y,X))).ExistAbstract(fp->cubePost_*fp->cubeOther_).SatCount(fp->preVars_.size()) << "\n";
            // std::cout << "states = " << term1.ExistAbstract(fp->cubePost_*fp->cubeOther_).SatCount(fp->preVars_.size()) << "\n";
            // std::cout << "states = " << term2.ExistAbstract(fp->cubePost_*fp->cubeOther_).SatCount(fp->preVars_.size()) << "\n";
//            sylvan::Bdd cubePre=sylvan::Bdd::VariablesCube(sort(fp->preVars_));
//            if (term2.ExistAbstract(cubePre)==sylvan::Bdd::bddOne()) {
//                std::cout << "Problem.";
//                if ((fp->apre(Y,X)).ExistAbstract(cubePre)==sylvan::Bdd::bddOne())
//                    std::cout << "Problem.";
//                if (term1.ExistAbstract(cubePre)==sylvan::Bdd::bddOne())
//                    std::cout << "Problem.";
//                if (right.ExistAbstract(cubePre)==sylvan::Bdd::bddOne())
//                    std::cout << "Problem.";
//                if (seqR.ExistAbstract(cubePre)==sylvan::Bdd::bddOne())
//                    std::cout << "Problem.";
//            }
            /* debug ends */
            /* add the recently added state-input pairs to the controller */
            N = term2 & (!(C.ExistAbstract(fp->cubePost_*fp->cubeOther_)));
            C |= N;
            if (remPairs.size()==0) {
                XX= term2;
                /* debug */
                // std::cout << "states = " << XX.ExistAbstract(fp->cubePost_*fp->cubeOther_).SatCount(fp->preVars_.size()) << "\n";
                /* debug ends */
            } else {
                const_arg_recursive_rabin arg_const_new = {
                    accl_on,
                    M, /* the bound on the iteration count for memorizing the BDDs from the past iterations */
                    depth+1,
                    remPairs,
                    initial_seed,
                    verbose

                };
                nconst_arg_recursive_rabin arg_nconst_new = {(seqR & nR).ExistAbstract(fp->cubePost_*fp->cubeOther_),
                    term2.ExistAbstract(fp->cubePost_*fp->cubeOther_),
                    indexRP,
                    indexY,
                    indexX,
                    hist_Y,
                    hist_X
                };
                XX= CALL(RabinRecurse, fp, &C, &arg_const_new, &arg_nconst_new);
            }
//             /* debug */
//             if (k>30) {
//                 std::cout << "Problem.";
// //                store_bdd_in_file((X), fp->preVars_, fp->postVars_, "X.txt");
// //                store_bdd_in_file((XX), fp->preVars_, fp->postVars_, "XX.txt");
// //                store_bdd_in_file(((!XX)&X), fp->preVars_, fp->postVars_, "XnXX.txt");
//             }
// //            if (!(X<=XX)) {
//             if ((X.ExistAbstract(fp->cubePost_*fp->cubeOther_) & (!XX.ExistAbstract(fp->cubePost_*fp->cubeOther_))) != (sylvan::Bdd::bddZero())) {
//                 sylvan::Bdd EX=X.ExistAbstract(fp->cubePost_*fp->cubeOther_);
//                 sylvan::Bdd EXX=XX.ExistAbstract(fp->cubePost_*fp->cubeOther_);
//             // if (!(X.ExistAbstract(fp->cubePost_*fp->cubeOther_)<=XX.ExistAbstract(fp->cubePost_*fp->cubeOther_))) {
//                 std::cout << "Non-monotonic mu fixpoint. Stop.";
//                 if (X<=XX)
//                     std::cout << "Weird\n\n";
//                 print_bdd_info(X.ExistAbstract(fp->cubePost_*fp->cubeOther_), fp->preVars_);
//                 print_bdd_info(XX.ExistAbstract(fp->cubePost_*fp->cubeOther_), fp->preVars_);
//                store_bdd_in_file((X), fp->preVars_, fp->postVars_, "X.txt");
//                store_bdd_in_file((XX), fp->preVars_, fp->postVars_, "XX.txt");
//                store_bdd_in_file(((!XX)&X), fp->preVars_, fp->postVars_, "XnXX.txt");
//                print_bdd_info(((!XX.ExistAbstract(fp->cubePost_*fp->cubeOther_))& X.ExistAbstract(fp->cubePost_*fp->cubeOther_)), fp->preVars_);
// //                exit(0);
//             }
// //            cubePre=sylvan::Bdd::VariablesCube(sort(fp->preVars_));
// //            if (X.ExistAbstract(cubePre)==sylvan::Bdd::bddOne())
// //                std::cout << "Problem.";
// //            if (k>(fp->nodes_).SatCount(fp->preVars_.size())) {
// //                std::cout << "Problem.";
// //            }
//             /* debug ends */
            if (accl_on)
                indexX->pop_back();
        }
        YY=XX;
        /* debug */
//        if (!(YY<=Y)) {
//            std::cout << "Non-monotonic nu fixpoint. Stop.";
//            exit(0);
//        }
//        if (j>(fp->nodes_).SatCount(fp->preVars_.size())) {
//            std::cout << "Problem.";
//        }
        /* debug ends */
        if (accl_on) {
            if (accl_on && fp->check_threshold(*indexY,M-1) && fp->check_threshold(*indexX,M-1)) {
//                if ((*hist_X)
//                    [fp->lexi_order(*indexRP,fp->RabinPairs_.size()-1)]
//                    [fp->to_dec(M,fp->pad_zeros(*indexY,remPairs.size()))]
//                    [depth-1]
//                    [std::min((*indexX)[0]+1,M-1)] <= XX.ExistAbstract(fp->cubePost_*fp->cubeOther_)*fp->tr_) {
//                    (*hist_X)
//                    [fp->lexi_order(*indexRP,fp->RabinPairs_.size()-1)]
//                    [fp->to_dec(M,fp->pad_zeros(*indexY,remPairs.size()))]
//                    [depth-1]
//                    [std::min((*indexX)[0]+1,M-1)]
//                    =XX.ExistAbstract(fp->cubePost_*fp->cubeOther_)*fp->tr_;
//                }
//                (*hist_X)[fp->lexi_order(*indexRP,fp->RabinPairs_.size()-1)][fp->to_dec(M,fp->pad_zeros(*indexY,remPairs.size()))][depth-1]=XX.ExistAbstract(fp->cubePost_*fp->cubeOther_)*fp->tr_;
                if ((*hist_X)[depth-1]\
                    [fp->findRank(fp->RabinPairs_.size(),*indexRP)]\
                    [std::min((*indexX)[0]+1,M-1)]\
                    [fp->to_dec(M,*indexY)] <= (XX.ExistAbstract(fp->cubePost_*fp->cubeOther_)*fp->tr_)) {
                    (*hist_X)[depth-1]\
                        [fp->findRank(fp->RabinPairs_.size(),*indexRP)]\
                        [std::min((*indexX)[0]+1,M-1)]\
                    [fp->to_dec(M,*indexY)] = (XX.ExistAbstract(fp->cubePost_*fp->cubeOther_)*fp->tr_);
                }
            }
            indexY->pop_back();
        }
    }
    *controller = C;
    if (accl_on) {
        if (accl_on && fp->check_threshold(*indexY,M-1) && fp->check_threshold(*indexX,M-1)) {
//            if (YY.ExistAbstract(fp->cubePost_*fp->cubeOther_)*fp->tr_ <=
//                (*hist_Y)
//                [fp->lexi_order(*indexRP,fp->RabinPairs_.size()-1)]
//                [fp->to_dec(M,fp->pad_zeros(*indexX,remPairs.size()+1))]
//                [depth-1]
//                [std::min((*indexY)[0]+1,M-1)]) {
//                (*hist_Y)
//                [fp->lexi_order(*indexRP,fp->RabinPairs_.size()-1)]
//                [fp->to_dec(M,fp->pad_zeros(*indexX,remPairs.size()+1))]
//                [depth-1]
//                [std::min((*indexY)[0]+1,M-1)]
//                =YY.ExistAbstract(fp->cubePost_*fp->cubeOther_)*fp->tr_;
//            }
//            (*hist_Y)[fp->lexi_order(*indexRP,fp->RabinPairs_.size()-1)][fp->to_dec(M,fp->pad_zeros(*indexX,remPairs.size()+1))][depth-1]=YY.ExistAbstract(fp->cubePost_*fp->cubeOther_)*fp->tr_;
            if ((YY.ExistAbstract(fp->cubePost_*fp->cubeOther_)*fp->tr_) <= (*hist_Y)[depth-1]\
                [fp->findRank(fp->RabinPairs_.size(),*indexRP)]\
                [std::min((*indexY)[0]+1,M-1)]\
                [fp->to_dec(M,*indexX)]) {
                (*hist_Y)[depth-1]\
                    [fp->findRank(fp->RabinPairs_.size(),*indexRP)]\
                    [std::min((*indexY)[0]+1,M-1)]\
                [fp->to_dec(M,*indexX)] = (YY.ExistAbstract(fp->cubePost_*fp->cubeOther_)*fp->tr_);
            }
        }
    }
    if (accl_on)
        indexRP->pop_back();

    YY=YY.ExistAbstract(fp->cubePost_*fp->cubeOther_);
    return YY.GetBDD();
}
} /* close namespace */
