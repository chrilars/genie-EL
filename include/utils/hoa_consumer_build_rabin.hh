//==============================================================================
//
// Author: Kaushik Mallik
//
// An interface to the cpphoafparser library for parsing a Rabin automaton in
// HOA format
//==============================================================================

#pragma once

#include <iostream>
#include <limits>
#include <map>
#include <stack>

#include "cpphoafparser/consumer/hoa_consumer.hh" // todo fix include
#include "cpphoafparser/parser/hoa_parser.hh" // todo fix include
#include <cstring> // todo fix include

#define UNUSED(x) (void)(x)

namespace cpphoafparser {

    struct transition {
        size_t state_id;
        HOAConsumer::label_expr::ptr label;
        size_t post_id;
    };

    struct rabin_data {
        size_t States;
        size_t Start;
        std::map<size_t, std::string> ap_id_map;
        std::vector<transition> Transitions;
        std::vector<std::vector<size_t>> acc_signature;
        std::vector<std::array<size_t, 2>> acc_pairs;
    };

    /**
 * A HOAConsumer implementation that works as an interface to convert a rabin automaton from HOA format
 * to the format used in Mascot-SDS.
 *
 */

    class HOAConsumerBuildRabin : public HOAConsumer {
    public:
        rabin_data *data_;

        /** Constructor, providing a reference to the output stream */
        HOAConsumerBuildRabin(rabin_data *data) : data_(data), out(std::cout) {}

        bool parserResolvesAliases() override {
            return false;
        }

        void notifyHeaderStart(const std::string &version) override {
            UNUSED(version);
            UNUSED(out);
        }

        void setNumberOfStates(unsigned int numberOfStates) override {
            data_->States = numberOfStates;
        }

        void addStartStates(const int_list &stateConjunction) override {
            if (stateConjunction.empty()) {
                throw std::runtime_error("no initial state specified");
            } else if (stateConjunction.size() > 1) {
                throw std::runtime_error("multiple initial states not allowed");
            }
            data_->Start = stateConjunction[0];
        }

        void addAlias(const std::string &name, label_expr::ptr labelExpr) override {
            UNUSED(name);
            UNUSED(labelExpr);
        }

        void setAPs(const std::vector<std::string> &aps) override {
            data_->ap_id_map.clear();
            for (size_t i = 0; i < aps.size(); i++) {
                data_->ap_id_map.insert({i, aps[i]});
            }
        }

        void setAcceptanceCondition(unsigned int numberOfSets, acceptance_expr::ptr accExpr) override {
            for (size_t i = 0; i < numberOfSets; i++) {
                std::vector<size_t> acc_label;
                data_->acc_signature.push_back(acc_label);
            }
            std::stack<acceptance_expr::ptr> nodes;
            nodes.push(accExpr);
            while (nodes.empty()) {
                acceptance_expr::ptr curr_node = nodes.top();
                nodes.pop();
                if (curr_node->isOR()) {
                    nodes.push(curr_node->getLeft());
                    nodes.push(curr_node->getRight());
                } else if (curr_node->isAND()) {
                    std::array<size_t, 2> pair;
                    acceptance_expr::ptr lchild = curr_node->getLeft();
                    acceptance_expr::ptr rchild = curr_node->getRight();
                    if (!lchild->isAtom() || !rchild->isAtom()) {
                        throw std::runtime_error("not a valid rabin specification.");
                    } else {
                        if (lchild->getAtom().getType() == AtomAcceptance::TEMPORAL_FIN) {
                            pair[1] = lchild->getAtom().getAcceptanceSet();
                            if (rchild->getAtom().getType() != AtomAcceptance::TEMPORAL_INF) {
                                throw std::runtime_error("not a valid rabin specification.");
                            } else {
                                pair[0] = rchild->getAtom().getAcceptanceSet();
                            }
                        } else if (lchild->getAtom().getType() == AtomAcceptance::TEMPORAL_INF) {
                            pair[0] = lchild->getAtom().getAcceptanceSet();
                            if (rchild->getAtom().getType() != AtomAcceptance::TEMPORAL_FIN) {
                                throw std::runtime_error("not a valid rabin specification.");
                            } else {
                                pair[1] = rchild->getAtom().getAcceptanceSet();
                            }
                        } else {
                            throw std::runtime_error("not a valid rabin specification.");
                        }
                    }
                    data_->acc_pairs.push_back(pair);
                }
            }
        }

        void
        provideAcceptanceName(const std::string &name, const std::vector<IntOrString> &extraInfo) override {
            if (std::strcmp(name.c_str(), "Rabin") != 0) {
                throw std::runtime_error("the input automaton is not rabin automaton");
            }
            UNUSED(extraInfo);
        }

        void setName(const std::string &name) override {
            UNUSED(name);
        }

        void setTool(const std::string &name, std::shared_ptr<std::string> version) override {
            UNUSED(name);
            UNUSED(version);
        }

        void addProperties(const std::vector<std::string> &properties) override {
            UNUSED(properties);
        }

        void addMiscHeader(const std::string &name, const std::vector<IntOrString> &content) override {
            UNUSED(name);
            UNUSED(content);
        }

        void notifyBodyStart() override {}

        void addState(unsigned int id,
                      std::shared_ptr<std::string> info,
                      label_expr::ptr labelExpr,
                      std::shared_ptr<int_list> accSignature) override {
            if (accSignature) {
                for (unsigned int acc: *accSignature) {
                    data_->acc_signature[acc].push_back(id);
                }
            }
            UNUSED(info);
            UNUSED(labelExpr);
        }

        void addEdgeImplicit(unsigned int stateId,
                             const int_list &conjSuccessors,
                             std::shared_ptr<int_list> accSignature) override {
            UNUSED(stateId);
            UNUSED(conjSuccessors);
            UNUSED(accSignature);
        }

        void addEdgeWithLabel(unsigned int stateId,
                              label_expr::ptr labelExpr,
                              const int_list &conjSuccessors,
                              std::shared_ptr<int_list> accSignature) override {
            if (conjSuccessors.size() != 1) {
                throw std::runtime_error("the rabin automaton must be deterministic");
            }
            transition t;
            t.state_id = stateId;
            t.label = labelExpr;
            t.post_id = conjSuccessors[0];

            data_->Transitions.push_back(t);

            UNUSED(accSignature);
        }

        void notifyEndOfState(unsigned int stateId) override {
            UNUSED(stateId);
        }

        void notifyEnd() override {}

        void notifyAbort() override {}

        void notifyWarning(const std::string &warning) override {
            std::cerr << "Warning: " << warning << std::endl;
        }

    private:
        /** Reference to the output stream */
        std::ostream &out;
    };

}// namespace cpphoafparser
