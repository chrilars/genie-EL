/*
* TransitionSystem.hh
*
*/

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <stdexcept>
#include <sylvan_obj.hpp>
#include <vector>

namespace fairsyn {
    struct thread {
        std::vector<uint32_t> state_vars;
        std::vector<uint32_t> state_vars_post;
        // std::vector<uint32_t> action_vars;
        sylvan::Bdd states;
        sylvan::Bdd transitions;
    };
    /*
    * class: Arena
    *
    *
    * a representation of the finite two-player turn-based game arena with edge fairness condition
    *
    */
    class Arena {
    public:
        /* var: nodes_
        * stores the nodes in a BDD */
        sylvan::Bdd nodes_;
        /* var: sys_nodes_
         * stores the system nodes in a BDD */
        sylvan::Bdd sys_nodes_;
        /* var: env_nodes_
        * stores the environment nodes in a BDD */
        sylvan::Bdd env_nodes_;
        //    /* var: inputSpace_
        //     * stores the input space BDD */
        //    sylvan::Bdd inputSpace_;
        /* var: preVars_
        * stores the "pre" node variable indices */
        std::vector<uint32_t> preVars_;
        //    /* var: inputVars_
        //     * stores the input variable indices */
        //    std::vector<uint32_t> inputVars_;
        /* var: postVars_
      * stores the "post" node variable indices */
        std::vector<uint32_t> postVars_;
        /* var: tr_
       * Transition BDD */
        sylvan::Bdd tr_;
        /* var: live_
          * the live transitions (subset of the transition relation) */
        sylvan::Bdd live_;


    public:
        /*
    * constructor: Arena
    *
    *  construct Arena from a given product of thread interfaces
    */
        Arena(std::set<uint32_t> &all_variables,
              std::vector<thread *> &thread_list,
              std::vector<uint32_t> &resource_vars,
              sylvan::Bdd resource_states,
              std::vector<uint32_t> &resource_vars_post,
              std::vector<uint32_t> &input_action_vars,
              std::vector<uint32_t> &output_action_vars,
              sylvan::Bdd input_action_states,
              sylvan::Bdd output_action_states,
              std::map<std::string, sylvan::Bdd> &progress_pred_to_bdd) {
            for (size_t i = 0; i < thread_list.size(); i++) {
                for (auto j = thread_list[i]->state_vars.begin(); j != thread_list[i]->state_vars.end(); ++j) {
                    all_variables.insert(*j);
                    preVars_.push_back(*j);
                }
                for (auto j = thread_list[i]->state_vars_post.begin(); j != thread_list[i]->state_vars_post.end(); ++j) {
                    all_variables.insert(*j);
                    postVars_.push_back(*j);
                }
            }
            for (auto j = resource_vars.begin(); j != resource_vars.end(); ++j) {
                all_variables.insert(*j);
                preVars_.push_back(*j);
            }
            for (auto j = resource_vars_post.begin(); j != resource_vars_post.end(); ++j) {
                all_variables.insert(*j);
                postVars_.push_back(*j);
            }
            /* NOTE: what to do with the action_vars: they go to preVars_ or postVars_? Or do they remain separate? I think they should remain separate (implicitly appears on the transitions). */
            for (auto j = input_action_vars.begin(); j != input_action_vars.end(); ++j) {
                all_variables.insert(*j);
            }
            for (auto j = output_action_vars.begin(); j != output_action_vars.end(); ++j) {
                all_variables.insert(*j);
            }
            /* ----------- the state space -----------------
        * state space of the product is the product of the state spaces of the threads and the resource state space */
            nodes_ = sylvan::Bdd::bddOne();
            for (size_t i = 0; i < thread_list.size(); i++) {
                nodes_ &= thread_list[i]->states;
            }
            nodes_ &= resource_states;
            // print_bdd_info(nodes_, thread_list[0]->state_vars);
            /* ------------ the transitions -----------------*/
            std::vector<uint32_t> thread_vars, thread_vars_post; /* bdd variables for representing which of the threads moves in a transition
       When all the variables are 0, it means that no thread moves (for example, from the system nodes). */
            for (uint32_t i = 1; i <= std::ceil(std::log2(thread_list.size() + 1)); i++) {
                uint32_t new_var_id = *all_variables.rbegin() + 1;
                thread_vars.push_back(new_var_id);
                preVars_.push_back(new_var_id);
                all_variables.insert(new_var_id);
            }
            for (uint32_t i = 1; i <= std::ceil(std::log2(thread_list.size() + 1)); i++) {
                uint32_t new_var_id = *all_variables.rbegin() + 1;
                thread_vars_post.push_back(new_var_id);
                postVars_.push_back(new_var_id);
                all_variables.insert(new_var_id);
            }
            sylvan::BddSet thread_vars_bdd = sylvan::BddSet::fromVector(thread_vars);
            sylvan::BddSet thread_vars_post_bdd = sylvan::BddSet::fromVector(thread_vars_post);
            /* build the product transition with the following constraints:
            * - at a time, only one thread moves to the next step and the rest remains standstill
            * - the index of the moving thread is recorded in the thread_vars_post bdd variables
            */
            tr_ = sylvan::Bdd::bddZero();
            for (size_t i = 0; i < thread_list.size(); i++) {
                /* the bdd representing the transitions where only thread i moves */
                sylvan::Bdd T = sylvan::Bdd::bddOne();
                /* thread i moves */
                T &= elementToBdd(i + 1, thread_vars_post_bdd, thread_vars_post.size()); /* select thread i; the "+1" is to account for the fact that the decimal value "0" is reserved for the case when no thread moves. */
                T &= thread_list[i]->transitions;
                /* rest of the threads remain at standstill */
                for (size_t j = 0; j < thread_list.size(); j++) {

                    if (i == j) {
                        continue;
                    }

                    for (size_t k = 0; k < thread_list[j]->state_vars.size(); k++) {
                        /* b represents a bdd of the form (x <=> x'), where x and x' are the k-th pre and post state variables of the thread j */
                        sylvan::Bdd b = sylvan::Bdd::bddVar(thread_list[j]->state_vars[k]) & sylvan::Bdd::bddVar(thread_list[j]->state_vars_post[k]);
                        b |= (!sylvan::Bdd::bddVar(thread_list[j]->state_vars[k])) & (!sylvan::Bdd::bddVar(thread_list[j]->state_vars_post[k]));
                        T &= b;
                    }
                }

                tr_ |= T;
            }
            tr_ &= nodes_;
            /* debug */
            // std::vector<uint32_t> actionVars=input_action_vars;
            // actionVars.push_back(output_action_vars[0]);
            // print_bdd_info(tr_, preVars_, actionVars, postVars_);
            // print_bdd_info(thread_list[0]->transitions, thread_list[0]->state_vars, actionVars, thread_list[0]->state_vars_post);
            // print_bdd_info(tr_, sort(preVars_), sort(postVars_));
            /* debug ends */
            // tr_ &= nodes_;
            // print_bdd_info(tr_, preVars_, postVars_);
            /* return the mapping from progress_i symbols to the state predicates recording the progress in the individual threads */
            progress_pred_to_bdd.clear();
            for (size_t i = 0; i < thread_list.size(); i++) {
                std::string s = "progress_";
                s += std::to_string(i);
                // progress_pred_to_bdd.insert({s, nodes_ & elementToBdd(i+1, thread_vars_bdd, thread_vars.size())}); /* mark progress in thread i; the "+1" is to account for the fact that the decimal value "0" is reserved for the case when no thread moves. */
                progress_pred_to_bdd.insert({s, nodes_.Permute(preVars_, postVars_) & elementToBdd(i + 1, thread_vars_post_bdd, thread_vars_post.size())}); /* mark progress in thread i; the "+1" is to account for the fact that the decimal value "0" is reserved for the case when no thread moves. */
            }
            /* debug */
            // print_bdd_info(tr_&progress_pred_to_bdd["progress_0"], preVars_, postVars_);
            // std::vector<uint32_t> temp=preVars_;
            // temp.push_back(output_action_vars[0]);
            // print_bdd_info(tr_, temp, postVars_);
            /* debug ends */
            /* ---------------- the game (splitting states into system and environment states) ------------ */
            uint32_t system_nodes = *all_variables.rbegin() + 1; /* bdd variable that is 1 for system nodes and 0 for the environment nodes */
            sylvan::Bdd system_nodes_bdd = sylvan::Bdd::bddVar(system_nodes);
            sylvan::Bdd environment_nodes_bdd = !system_nodes_bdd;
            sys_nodes_ = nodes_ & system_nodes_bdd;
            env_nodes_ = nodes_ & (!system_nodes_bdd);
            preVars_.push_back(system_nodes);
            all_variables.insert(system_nodes);
            uint32_t system_nodes_post = *all_variables.rbegin() + 1; /* the post bdd variable */
            sylvan::Bdd system_nodes_post_bdd = sylvan::Bdd::bddVar(system_nodes_post);
            sylvan::Bdd environment_nodes_post_bdd = !system_nodes_post_bdd;
            postVars_.push_back(system_nodes_post);
            all_variables.insert(system_nodes_post);
            /* encode the interplay between the system and the environment */
            sylvan::Bdd game_rule = sylvan::Bdd::bddZero();
            /* from the system nodes, the game can move to one of these:
            * (1) an input successor
            * (2) the environment node which has the same valuations for the rest of the bdd variables */
            /* rule (1) */
            game_rule |= (system_nodes_bdd & input_action_states & (!output_action_states) & system_nodes_post_bdd);

            /* debug */
            // print_bdd_info(game_rule, preVars_, actionVars, postVars_);
            /* debug ends */

            /* rule (2) */
            sylvan::Bdd b = sylvan::Bdd::bddOne();
            b = system_nodes_bdd & environment_nodes_post_bdd;
            for (size_t k = 0; k < preVars_.size(); k++) {
                if (preVars_[k] == system_nodes) {
                    continue;
                }
                /* b represents a bdd of the form (x <=> x'), where x and x' are the k-th pre and post state variables of the product state space */
                sylvan::Bdd c = sylvan::Bdd::bddVar(preVars_[k]) & sylvan::Bdd::bddVar(postVars_[k]);
                c |= (!sylvan::Bdd::bddVar(preVars_[k])) & (!sylvan::Bdd::bddVar(postVars_[k]));
                b &= c;
            }
            game_rule |= b;

            /* debug */
            // print_bdd_info(game_rule, preVars_, actionVars, postVars_);
            /* debug ends */

            /* from the envrionment nodes, the game can move along one of the output successors, and the next state is a sysetm node */
            game_rule |= (environment_nodes_bdd & output_action_states & (!input_action_states) & system_nodes_post_bdd);

            /* debug */
            // temp=preVars_;
            // temp.push_back(output_action_vars[0]);
            // print_bdd_info(game_rule, temp, postVars_);
            // print_bdd_info(tr_, temp, postVars_);
            /* debug ends */

            /* debug */
            // print_bdd_info(game_rule, preVars_, actionVars, postVars_);
            // print_bdd_info(tr_, preVars_, actionVars, postVars_);
            /* debug ends */

            /* any system node can transition to an environment node without changing the other components */
            b = sylvan::Bdd::bddOne();
            b &= sylvan::Bdd::bddVar(system_nodes) & (!sylvan::Bdd::bddVar(system_nodes_post));
            for (size_t k = 0; k < preVars_.size(); k++) {
                if (preVars_[k] == system_nodes) {
                    continue;
                }
                /* b represents a bdd of the form (x XOR x'), where x and x' are the k-th pre and post state variables of the product state space */
                sylvan::Bdd c = sylvan::Bdd::bddVar(preVars_[k]) & sylvan::Bdd::bddVar(postVars_[k]);
                c |= (!sylvan::Bdd::bddVar(preVars_[k])) & (!sylvan::Bdd::bddVar(postVars_[k]));
                b &= c;
            }
            tr_ |= b;
            tr_ &= nodes_;

            /* finally, incorporate the system-environment interaction in the transitions */
            tr_ &= game_rule;
            /* debug */
            // print_bdd_info(tr_, preVars_, actionVars, postVars_);
            /* debug ends */
            // print_bdd_info(game_rule, preVars_, postVars_);
            // store_bdd_in_file(tr_, preVars_, postVars_, "tr.txt");
            /* the live transitions are all the environment transitions */
            live_ = tr_ & environment_nodes_bdd;
            /* sort the list of variables */
            std::sort(preVars_.begin(), preVars_.end());
            std::sort(postVars_.begin(), postVars_.end());
            // /* debug */
            // print_bdd_info(tr_, preVars_, postVars_);
            // print_bdd_info(progress_pred_to_bdd["progress_0"], preVars_);
            // /* debug ends */
        }
        /*
    * constructor: Arena
    */
        Arena(std::set<uint32_t> &all_variables,
              std::vector<size_t> &nodes,
              std::vector<size_t> &sys_nodes,
              std::vector<size_t> &env_nodes,
              std::vector<std::vector<size_t>> &transitions,
              std::vector<std::vector<size_t>> &live_edges) {
            /* number of bdd variables required for the nodes */
            uint32_t nvars = std::ceil(std::log2(nodes.size()));
            if (all_variables.size() == 0) {
                /* the "pre" bdd variables */
                for (uint32_t i = 0; i < nvars; i++) {
                    uint32_t new_var_id = i;
                    preVars_.push_back(new_var_id);
                    all_variables.insert(new_var_id);
                }
            } else {
                /* the "pre" bdd variables */
                uint32_t max_var_id = *(all_variables.rbegin());
                for (uint32_t i = 1; i <= nvars; i++) {
                    uint32_t new_var_id = max_var_id + i;
                    preVars_.push_back(new_var_id);
                    all_variables.insert(new_var_id);
                }
            }
            uint32_t max_var_id = *(all_variables.rbegin());

            /* the "post" bdd variables */
            for (uint32_t i = 1; i <= nvars; i++) {
                uint32_t new_var_id = max_var_id + i;
                postVars_.push_back(new_var_id);
                all_variables.insert(new_var_id);
            }
            /* the bdd representing the set of nodes
        * given in terms of the "pre" variables */
            sylvan::BddSet preBddVars = sylvan::BddSet::fromVector(preVars_);
            nodes_ = setToBdd(nodes, preBddVars, nvars);
            sys_nodes_ = setToBdd(sys_nodes, preBddVars, nvars);
            env_nodes_ = setToBdd(env_nodes, preBddVars, nvars);
            /* the bdd representing the transition relation */
            tr_ = sylvan::Bdd::bddZero();
            sylvan::BddSet postBddVars = sylvan::BddSet::fromVector(postVars_);
            for (size_t i = 0; i < transitions.size(); i++) {
                /* the i-th element of the transitions vector is the vector of successors from the i-th vertex */
                sylvan::Bdd cur_node = elementToBdd(i, preBddVars, nvars);
                sylvan::Bdd succ_node = setToBdd(transitions[i], postBddVars, nvars);
                tr_ += (cur_node & succ_node);
            }

            /*for asserting whether BDDs formed are right or not */
            for (size_t i = 0; i < transitions.size(); i++) {
                sylvan::Bdd cur_node = elementToBdd(i, preBddVars, nvars);

                for (size_t j = 0; j < transitions[i].size(); j++) {
                    sylvan::Bdd temp = sylvan::Bdd::bddZero();
                    sylvan::Bdd succ_node = elementToBdd(transitions[i][j], postBddVars, nvars);
                    temp = (cur_node & succ_node);
                    sylvan::Bdd intermediate = temp & tr_;
                    if (temp != intermediate) {
                        std::cout << " BDDs not formed well \n";
                        exit(0);
                    }
                }
            }


            /* the bdd representing the live edges */
            live_ = sylvan::Bdd::bddZero();
            for (size_t i = 0; i < live_edges.size(); i++) {
                /* the i-th element of the live transitions vector is the vector of live edge successors from the i-th vertex */
                sylvan::Bdd cur_node = elementToBdd(i, preBddVars, nvars);
                sylvan::Bdd succ_node = setToBdd(live_edges[i], postBddVars, nvars);
                live_ += (cur_node & succ_node);
            }

            /* for asserting that the system nodes and the environment nodes form a partition only */
            if (((sys_nodes_ & env_nodes_) != sylvan::Bdd::bddZero()) |
                ((sys_nodes_ | env_nodes_) != nodes_)) {
                std::cout << "The system and the environment nodes do not form a partition. Exiting.\n";
                exit(0);
            }

            /* for asserting that the live edges are from env vertices only */
            if ((live_ & sys_nodes_) != sylvan::Bdd::bddZero()) {
                std::cout << "Some live edges start from the system nodes. Exiting.\n";
                exit(0);
            }

            /* for asserting that transitions are valid */
            if (!(tr_ <= (nodes_ & nodes_.Permute(preVars_, postVars_)))) {
                std::cout << "Some transitions are not valid. Exiting.\n";
                exit(0);
            }

            /* for asserting that the live edges are subsets of the transitions */
            if (!(live_ <= tr_)) {
                std::cout << "Some live edges are not valid transitions. Exiting.\n";
                exit(0);
            }
        }
        /*
    * Default constructor
    */
        Arena() {
            nodes_ = sylvan::Bdd::bddZero();
            sys_nodes_ = sylvan::Bdd::bddZero();
            env_nodes_ = sylvan::Bdd::bddZero();
            tr_ = sylvan::Bdd::bddZero();
            live_ = sylvan::Bdd::bddZero();
        }
        /* function: sort */
        template <class T>
        inline std::vector<T> sort(const std::vector<T> &vec) {
            std::vector<T> sorted = vec;
            std::sort(sorted.begin(), sorted.end());
            return sorted;
        }

    private:
        template<typename T>
        sylvan::Bdd setToBdd(const T& set,
                             const sylvan::BddSet& bddVars,
                             const size_t nvars) {
            sylvan::Bdd set_bdd=sylvan::Bdd::bddZero();
            using std::begin;
            using std::end;
            for (auto i=begin(set); i!=end(set); ++i) {
                /* the current element */
                uint32_t x=*i;
                /* compute the binary representation of the current element */
                std::vector<uint8_t> bin(nvars, 0);
                size_t p=0;
                for (; x; x/=2, p++) {
                    bin[p]=0+x%2;
                }
                /* update the set bdd */
                set_bdd=set_bdd.UnionCube(bddVars, bin);
            }
            return set_bdd;
        }
        
        sylvan::Bdd elementToBdd(const size_t& x,
                                 const sylvan::BddSet& bddVars,
                                 const size_t nvars) {
            size_t y=x;
            /* debug */
            std::vector<uint32_t> v = bddVars.toVector();
            /* debug ends */
            /* compute the binary representation of the current element */
            std::vector<uint8_t> bin(nvars, 0);
            size_t p=0;
            for (; y; y/=2, p++) bin[p]=0+y%2;
            /* return the cube */
            return sylvan::Bdd::bddCube(bddVars, bin);
        }
        
    }; /* close class definition */
}// namespace fairsyn
