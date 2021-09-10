/*
 * RabinAutomaton.hh
 *
 *  created on: 28.01.2021
 *      author: kaushik mallik
 *
 */

#ifndef RABINAUTOMATON_HH_
#define RABINAUTOMATON_HH_

#include <math.h>
#include <stdlib.h>

#include "Helper.hh"
#include "MascotHelper.hh"
#include "SymbolicSet.hh"
#include "utils/hoa_consumer_build_rabin.hh"

/* class: RabinAutomaton
 *
 * symbolic (bdd) implementation of a Rabin automaton
 *
 */
namespace mascot {
    template <class UBDD>
    class RabinAutomaton {
    public:
        /* var: base_
         *  the bdd manager */
        UBDD base_;
        /* var: stateSpace_ */
        SymbolicSet<UBDD> *stateSpace_;
        /* var: inputSpace_ */
        SymbolicSet<UBDD> *inputSpace_;
        /* var: stateSpacePost_ */
        SymbolicSet<UBDD> *stateSpacePost_;
        /* var: transitions_
         *  the bdd representing the transition relation as sets of tuples (s,x',s'), where s is pre-state of the rabin automaton, x' is post-state of the symbolic model, and s' is the post-state of the rabin automaton */
        UBDD transitions_;
        /* var: numRabinPairs_
         *  the number of rabin pairs */
        size_t numRabinPairs_;
        /* var: RabinPairs_
         * UBDD vector[numRabinPairs_][2] containing the rabin pairs {(G_i,~R_i)}, where the Rabin condition is given as: \/_i ([]<>G_i & <>[]~R_i)
         * and each G_i,~R_i are given a bdd representing a set of states of the automaton */
        std::vector<rabin_pair_<UBDD>> RabinPairs_;

    public:
        /* constructor: RabinAutomaton
         * constructs the bdd representation of the Rabin automaton from a given list representation of the states and the transitions
         *
         * Input:
         * base    - the UBDD manager
         * ns       - number of states
         * inputSpace - symbolic set representing the input space
         * inputs   - vector containing the inputs; index should match with the transitions (predicate over the post states of the symbolic model)
         * transitions - vector containing the (directed) transitions (as arrays of size 3, containing the source vertex, the input, and the end vertex)
         * RabinPairs - vector containing the list of Rabin pairs (as pairs of Gi,Ri, where each Gi and Ri are vectors of state indices) */
        RabinAutomaton(UBDD &base,
                       const size_t ns,
                       const SymbolicSet<UBDD> &inputSpace,
                       const std::vector<UBDD> inputs,
                       const std::vector<std::array<size_t, 3>> transitions,
                       const std::vector<std::array<std::vector<size_t>, 2>> RabinPairs) {
            base_ = base;
            /* create a symbolic set representing the state space of the rabin automaton
             * the symbolic set is created by discretizing the range [1,ns] on the real line with eta=1 */
            size_t sdim = 1; /* we model the states as discrete points on the real line */
            std::vector<double> slb = {0};
            std::vector<double> sub = {static_cast<double>(ns - 1)};
            std::vector<double> seta = {1};
            int grid_alignment = 0; /* the origin is a grid point */
            stateSpace_ = new SymbolicSet<UBDD>(base_, sdim, slb, sub, seta, grid_alignment);
            stateSpace_->addGridPoints();
            inputSpace_ = new SymbolicSet<UBDD>(inputSpace);
            inputSpace_->addGridPoints();
            stateSpacePost_ = new SymbolicSet<UBDD>(*stateSpace_, 1);
            stateSpacePost_->addGridPoints();
            /* build the transitions bdd */
            transitions_ = base_.zero();
            for (size_t i = 0; i < transitions.size(); i++) {
                UBDD cube = base_.one();
                std::vector<size_t> v;
                v = {transitions[i][0]};
                cube &= stateSpace_->elementToMinterm(v);
                cube &= inputs[transitions[i][1]];
                v = {transitions[i][2]};
                cube &= stateSpacePost_->elementToMinterm(v);
                transitions_ |= cube;
            }
            /* build the rabin pairs in terms of the pre variables */
            numRabinPairs_ = RabinPairs.size();
            for (size_t i = 0; i < RabinPairs.size(); i++) {
                rabin_pair_<UBDD> pair;
                pair.rabin_index_ = i;
                pair.G_ = base_.zero();
                pair.nR_ = base_.zero();
                std::vector<size_t> v;
                /* first, create theUBDD for the G sets */
                for (size_t j = 0; j < RabinPairs[i][0].size(); j++) {
                    v = {RabinPairs[i][0][j]};
                    pair.G_ |= stateSpace_->elementToMinterm(v);
                }
                /* second, create theUBDD for the COMPLEMENT OF the R sets */
                for (size_t j = static_cast<size_t>(stateSpace_->getFirstGridPoint()[0]); j <= static_cast<size_t>(stateSpace_->getLastGridPoint()[0]); j++) {
                    bool isInR = false;
                    for (size_t k = 0; k < RabinPairs[i][1].size(); k++) {
                        if (RabinPairs[i][1][k] == j) {
                            isInR = true;
                            break;
                        }
                    }
                    if (!isInR) {
                        v = {j};
                        pair.nR_ |= stateSpace_->elementToMinterm(v);
                    }
                }
                RabinPairs_.push_back(pair);
            }
        }
        /* constructor: RabinAutomaton
         * constructs the bdd representation of the Rabin automaton from a given list representation of the states and the transitions
         *
         * Input:
         * base    - the UBDD manager
         * inputSpace - symbolic set representing the input space
         * hoaf     - the file where the rabin automaton has been given in HOA format
        */
        RabinAutomaton(UBDD &base,
                       const SymbolicSet<UBDD> &inputSpace,
                       std::map<std::string, UBDD> inputs_name_to_ubdd,
                       const std::string& filename) {
            base_ = base;
            std::ifstream file(filename);
            cpphoafparser::HOAConsumer::ptr consumer;
            cpphoafparser::rabin_data data;
            consumer.reset(new cpphoafparser::HOAConsumerBuildRabin(&data));
            cpphoafparser::HOAParser::parse(file, consumer);
            std::map<std::string, UBDD> inputs_id_to_ubdd;
            for (auto it = data.ap_id_map.begin(); it != data.ap_id_map.end(); ++it) {
                size_t key = it->first;
                inputs_id_to_ubdd.insert({std::to_string(key), inputs_name_to_ubdd[data.ap_id_map[key]]});
            }
            /* create a symbolic set representing the state space of the rabin automaton
             * the symbolic set is created by discretizing the range [1,ns] on the real line with eta=1 */
            size_t sdim = 1; /* we model the states as discrete points on the real line */
            std::vector<double> slb = {0};
            std::vector<double> sub = {static_cast<double>(data.States - 1)};
            std::vector<double> seta = {1};
            int grid_alignment = 0; /* the origin is a grid point */
            stateSpace_ = new SymbolicSet<UBDD>(base_, sdim, slb, sub, seta, grid_alignment);
            stateSpace_->addGridPoints();
            inputSpace_ = new SymbolicSet<UBDD>(inputSpace);
            inputSpace_->addGridPoints();
            stateSpacePost_ = new SymbolicSet<UBDD>(*stateSpace_, 1);
            stateSpacePost_->addGridPoints();
            /* build the transitions bdd */
            transitions_ = base_.zero();
            for (size_t i = 0; i < data.Transitions.size(); i++) {
                UBDD cube = base_.one();
                std::vector<size_t> v;
                v = {data.Transitions[i].state_id};
                cube &= stateSpace_->elementToMinterm(v);
                std::stack<cpphoafparser::HOAConsumer::label_expr::ptr> nodes;
                nodes.push(data.Transitions[i].label);
                //            cpphoafparser::HOAConsumer::label_expr::ptr curr_node = nodes.top();
                //            nodes.pop();
                cube &= label_to_ubdd(nodes, inputs_id_to_ubdd);
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
                v = {data.Transitions[i].post_id};
                cube &= stateSpacePost_->elementToMinterm(v);
                transitions_ |= cube;
            }
            /* build the rabin pairs in terms of the pre variables */
            numRabinPairs_ = data.acc_pairs.size();
            for (size_t i = 0; i < numRabinPairs_; i++) {
                rabin_pair_<UBDD> pair;
                pair.rabin_index_ = i;
                pair.G_ = base_.zero();
                pair.nR_ = base_.zero();
                std::vector<size_t> v;
                /* first, create theUBDD for the G sets */
                v = {data.acc_signature[data.acc_pairs[i][0]]};
                for (size_t j = 0; j < v.size(); j++) {
                    std::vector<size_t> temp = {v[j]};
                    pair.G_ |= stateSpace_->elementToMinterm(temp);
                }
                v = {data.acc_signature[data.acc_pairs[i][1]]};
                /* second, create theUBDD for the COMPLEMENT OF the R sets */
                for (size_t j = static_cast<size_t>(stateSpace_->getFirstGridPoint()[0]); j <= static_cast<size_t>(stateSpace_->getLastGridPoint()[0]); j++) {
                    bool isInR = false;
                    for (size_t k = 0; k < v.size(); k++) {
                        if (v[k] == j) {
                            isInR = true;
                            break;
                        }
                    }
                    if (!isInR) {
                        std::vector<size_t> temp = {j};
                        pair.nR_ |= stateSpace_->elementToMinterm(temp);
                    }
                }
                RabinPairs_.push_back(pair);
            }
        }
        /* function: label_to_ubdd */
        UBDD label_to_ubdd(std::stack<cpphoafparser::HOAConsumer::label_expr::ptr> &nodes,
                           std::map<std::string, UBDD> &inputs_id_to_ubdd) {
            cpphoafparser::HOAConsumer::label_expr::ptr curr_node = nodes.top();
            nodes.pop();
            if (curr_node->isAND()) {
                nodes.push(curr_node->getLeft());
                UBDD L = label_to_ubdd(nodes, inputs_id_to_ubdd);
                //            nodes.pop();
                nodes.push(curr_node->getRight());
                UBDD R = label_to_ubdd(nodes, inputs_id_to_ubdd);
                return (L & R);
                //            nodes.pop();
            } else if (curr_node->isOR()) { /* need to implement this */
                nodes.push(curr_node->getLeft());
                UBDD L = label_to_ubdd(nodes, inputs_id_to_ubdd);
                //            nodes.pop();
                nodes.push(curr_node->getRight());
                UBDD R = label_to_ubdd(nodes, inputs_id_to_ubdd);
                return (L | R);
                //            nodes.pop();
            } else if (curr_node->isNOT()) {
                return (!(inputs_id_to_ubdd[curr_node->getLeft()->toString()])) & inputSpace_->getSymbolicSet();
            } else if (curr_node->isAtom()) {
                return inputs_id_to_ubdd[curr_node->toString()];
            } else if (curr_node->isTRUE()) {
                return inputSpace_->getSymbolicSet();
            } else if (curr_node->isFALSE()) {
                return base_.zero();
            } else {
                throw std::runtime_error("Unknown operator in transition labels of the rabin automaton.");
            }
        }
        /* function:  writeToFile
     * save the transitions and the rabin pairs */
        void writeToFile(const char *foldername) {
            helper::checkMakeDir(foldername);
            /* create SymbolicSet representing the post stateSpace_*/
            /* make a copy of the SymbolicSet of the preSet */
            SymbolicSet<UBDD> post(*stateSpacePost_);
            /* create SymbolicSet representing the stateSpace_ x inputSpace_ */
            SymbolicSet<UBDD> sis(*stateSpace_, *inputSpace_);
            /* create SymbolicSet representing the stateSpace_ x inputSpace_ x stateSpace_ */
            SymbolicSet<UBDD> tr(sis, post);
            /* fill SymbolicSet with transition relation */
            tr.setSymbolicSet(transitions_);
            /* write the SymbolicSet to the file */
            std::string pathname = foldername;
            pathname += "/transitions.bdd";
            char Char[25];
            size_t Length = pathname.copy(Char, pathname.length() + 1);
            Char[Length] = '\0';
            tr.writeToFile(Char);
            /* create SymbolicSet-s representing the rabin pairs */
            SymbolicSet<UBDD> spec(*stateSpace_);
            for (size_t i = 0; i < RabinPairs_.size(); i++) {
                /* write the G sets */
                spec.setSymbolicSet(RabinPairs_[i].G_);
                std::string pathname = foldername;
                pathname += "/G";
                pathname += std::to_string(i);
                pathname += ".bdd";
                char Char[20];
                size_t Length = pathname.copy(Char, pathname.length() + 1);
                Char[Length] = '\0';
                spec.writeToFile(Char);
                /* write the nR sets */
                spec.setSymbolicSet(RabinPairs_[i].nR_);
                pathname = foldername;
                pathname += "/nR";
                pathname += std::to_string(i);
                pathname += ".bdd";
                Length = pathname.copy(Char, pathname.length() + 1);
                Char[Length] = '\0';
                spec.writeToFile(Char);
            }
        }

    private:
        /* function: dec2bin
         * decimal to binary conversion
         *
         * Input:
         * n - a number in decimal
         * Output:
         * b - a vector representing the bits of the binary representation, where b[0] is LSB */
        std::vector<bool> dec2bin(const size_t n) {
            std::vector<bool> b;
            if (n == 0) {
                b.push_back(0);
                return b;
            }
            size_t m = n;
            while (m > 0) {
                div_t d = div(m, 2);
                b.push_back(d.rem);
                m = d.quot;
            }
            return b;
        }
    }; /* close class def */
}// namespace mascot

#endif /* RABINAUTOMATON_HH_ */
