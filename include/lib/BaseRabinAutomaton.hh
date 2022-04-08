/** @file BaseRabinAutomaton.hh
*
*  @date 8.04.2022 (in progress)
*  @author Mateusz Rychlicki
*/

#pragma once

#include "../utils/hoa_consumer_build_rabin.hh"// todo fix that!
#include <functional>
#include <math.h>
#include <stdlib.h>
#include <vector>

namespace fairsyn {


    /** structure containing a single rabin pair */
    template <class UBDD>
    struct rabin_pair_ {
        size_t rabin_index_; /**< each rabin pair gets a unique index */
        UBDD G_;             /**< the states which are to be visited infinitely often */
        UBDD nR_;            /**< the complement of the states which are to be visited finitely often */
    };


    template <class UBDD>
    class BaseRabinAutomaton {
    public:
        UBDD base_;                                 /**< the bdd manager */
        UBDD transitions_;                          /**< @brief the bdd representing the transition relation
                                                     * @details the bdd representing the transition relation as sets of tuples (s,x',s'),
                                                     * where s is pre-state of the rabin automaton, x' is post-state of the symbolic model,
                                                     * and s' is the post-state of the rabin automaton */
        size_t numRabinPairs_;                      /**< the number of rabin pairs */
        std::vector<rabin_pair_<UBDD>> RabinPairs_; /**< @brief vector of sets of states of the automaton
                                                     * @details BDD vector[numRabinPairs_][2] containing the rabin pairs {(G_i,~R_i)},
                                                     * where the Rabin condition is given as:
                                                     * \/_i ([]<>G_i & <>[]~R_i)
                                                     * and each G_i,~R_i are given a bdd representing a set of states of the automaton */

        BaseRabinAutomaton(UBDD base):base_(base) {}

        virtual UBDD element_to_ubdd(size_t id) = 0;
        virtual UBDD element_to_ubdd_post(size_t id) = 0;

        void build_transitions(cpphoafparser::rabin_data &data,
                               std::map<std::string, UBDD> &inputs_id_to_ubdd,
                               UBDD isTrueCase) {
            transitions_ = base_.zero();
            size_t v;
            for (size_t i = 0; i < data.Transitions.size(); i++) {
                UBDD cube = base_.one();
                v = data.Transitions[i].state_id;
                cube &= element_to_ubdd(v);
                std::stack<cpphoafparser::HOAConsumer::label_expr::ptr> nodes;
                nodes.push(data.Transitions[i].label);
                //            cpphoafparser::HOAConsumer::label_expr::ptr curr_node = nodes.top();
                //            nodes.pop();
                cube &= label_to_ubdd(nodes, inputs_id_to_ubdd, isTrueCase);
                //            while (nodes.size()!=0) {
                //                cpphoafparser::HOAConsumer::label_expr::ptr curr_node = nodes.top();
                //                nodes.pop();
                //                if (curr_node->isAND()) {
                //                    nodes.push(curr_node->getLeft());
                //                    nodes.push(curr_node->getRight());
                //                } else if (curr_node->isOR()) { /* need to implement this */
                //                    throw std::runtime_error("Or not allowed in transition labels of the rabin automaton.");
                //                } else if (curr_node->isNOT()) {
                //                    cube &= (!(inputs_id_to_ubdd[curr_node->getLeft()->toString()])) & inputSpace_->getSymbolicSet();
                //                } else if (curr_node->isAtom()) {
                //                    cube &= inputs_id_to_ubdd[curr_node->toString()];
                //                } else if (curr_node->isTRUE()) {
                //                    cube &= inputSpace_->getSymbolicSet();
                //                } else if (curr_node->isFALSE()) {
                //                    cube &= base_.zero();
                //                }
                //            }
                //                v = {data.Transitions[i].post_id};
                cube &= element_to_ubdd_post(data.Transitions[i].post_id);
                transitions_ |= cube;
            }
        }

        /* function: label_to_ubdd */
        UBDD label_to_ubdd(std::stack<cpphoafparser::HOAConsumer::label_expr::ptr> &nodes,
                           std::map<std::string, UBDD> &inputs_id_to_ubdd,
                           UBDD isTrueCase) {
            cpphoafparser::HOAConsumer::label_expr::ptr curr_node = nodes.top();
            nodes.pop();
            if (curr_node->isAND()) {
                nodes.push(curr_node->getLeft());
                UBDD L = label_to_ubdd(nodes, inputs_id_to_ubdd, isTrueCase);
                //            nodes.pop();
                nodes.push(curr_node->getRight());
                UBDD R = label_to_ubdd(nodes, inputs_id_to_ubdd, isTrueCase);
                return (L & R);
                //            nodes.pop();
            } else if (curr_node->isOR()) { /* need to implement this */
                nodes.push(curr_node->getLeft());
                UBDD L = label_to_ubdd(nodes, inputs_id_to_ubdd, isTrueCase);
                //            nodes.pop();
                nodes.push(curr_node->getRight());
                UBDD R = label_to_ubdd(nodes, inputs_id_to_ubdd, isTrueCase);
                return (L | R);
                //            nodes.pop();
            } else if (curr_node->isNOT()) {
                return (!(inputs_id_to_ubdd[curr_node->getLeft()->toString()])) & isTrueCase;
            } else if (curr_node->isAtom()) {
                return inputs_id_to_ubdd[curr_node->toString()];
            } else if (curr_node->isTRUE()) {
                return isTrueCase;
            } else if (curr_node->isFALSE()) {
                return base_.zero();
            } else {
                throw std::runtime_error("Unknown operator in transition labels of the rabin automaton.");
            }
        }


        void build_rabin_pairs(std::function<std::vector<size_t> (size_t, size_t)> get_array) {
            std::vector<size_t> v;
            for (size_t i = 0; i < numRabinPairs_; i++) {
                rabin_pair_<UBDD> pair;
                pair.rabin_index_ = i;
                pair.G_ = base_.zero();
                /* first, create theUBDD for the G sets */
                v = get_array(i, 0);
                for (size_t j = 0; j < v.size(); j++) {
                    pair.G_ |= element_to_ubdd(v[j]);
                }
                /* second, create the UBDD for the COMPLEMENT OF the R sets */
                pair.nR_ = complement_of_R(get_array(i, 1));
                RabinPairs_.push_back(pair);
            }
        }

        virtual UBDD complement_of_R(const std::vector<size_t> &array) = 0;
    };
}// namespace fairsyn