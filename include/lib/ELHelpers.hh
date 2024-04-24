#pragma once

#include <algorithm>
#include <cassert>
#include <iostream>
#include <cmath>
#include <iterator>
#include <vector>
#include "structures.hh"

namespace ELHelpers {

    inline bool proper_subset(const std::vector<bool>& l1, const std::vector<bool>& l2) {
        // Check if l1 subset of l2
        // forall i in l1: if l1 then also l2
        // but there must exist at least one index in l2 which is not in l1
        bool flag = false;
        for (size_t i = 0; i < l1.size(); ++i) {
            if (l2[i] && !l1[i])
                flag = true;
            if (l1[i] && !l2[i])
                return false;
        }
        return flag;
    }

    inline void print_label(ZielonkaNode *z) {
        for (size_t i = 0; i < z->label.size(); ++i){
            if (z->label[i]) {
                std::cout << static_cast<char>('a' + i) << ", ";
            } else {
                std::cout << "   ";
            }
        } //std::cout << '\n';
    }


    inline std::vector<bool> label_difference(const std::vector<bool>& t, const std::vector<bool>& s) {
        // compute difference t-s, where s,t are vectors of colors
        size_t size = t.size();
        std::vector<bool> difference(size, false);
        for (size_t i = 0; i < size; ++i)
            difference[i] = t[i] && !s[i];
        return difference;
    }
    //template <class UBDD>
    //std::vector<UBDD> vector_difference(const std::vector<UBDD>& v1, const std::vector<UBDD>& v2) {
    //    // set of BDDs difference (v1 - v2)
    //    std::vector<UBDD> W;
    //    for (const auto& v : v1) {
    //        if (std::find(v2.begin(), v2.end(), v) == v2.end())
    //            W.push_back(v);
    //    }
    //    return W;
    //}

    //template <class UBDD>
    //std::vector<UBDD> vector_intersection(const std::vector<UBDD>& v1, const std::vector<UBDD>& v2) {
    //    // set of BDDS intersection
    //    std::vector<UBDD> W;
    //    for (const auto& v : v1) {
    //        if (std::find(v2.begin(), v2.end(), v) != v2.end())
    //            W.push_back(v);
    //    }
    //    return W;
    //}

    //template <class UBDD>
    //std::vector<UBDD> vector_union(const std::vector<UBDD>& v1, const std::vector<UBDD>& v2) {
    //    // set of BDDs union
    //    std::vector<UBDD> W = v2; // copy elements of v2 into W
    //    for (const auto& v : v1) {
    //        if (std::find(v2.begin(), v2.end(), v) != v2.end())
    //            W.push_back(v);
    //    }
    //    return W;
    //}

    //template <class UBDD>
    //std::vector<Vertex<UBDD>> gamma_inverse(std::vector<Color> colors, Arena* arena, bool subset_flag) {
    //    // return vertices in graph for which gamma is a subset of given set of colors
    //    std::vector<Vertex<UBDD>> vertices;
    //    for (const auto& v : arena->graph.vertices) {
    //        if (( subset_flag &&  subset(v.gamma, colors)) ||
    //            (!subset_flag && !subset(v.gamma, colors)))
    //            vertices.push_back(v);
    //    }
    //    return vertices;
    //}


    //template <class UBDD>
    //std::vector<Vertex<UBDD>> CPre(const UBDD& X, Arena* arena) {
    //    // TODO, Controlled predecessor function
    //    // CPre(X) = {v ∈ V∃ | E(v) ∩ X = ∅} ∪ {v ∈ V∀ | E(v) ⊆ X}
    //    // c1 = {v ∈ V∃ | E(v) ∩ X = ∅}
    //    // c2 = {v ∈ V∀ | E(v) ⊆ X}

    //    // c1(X): vertices v in existential vertices,
    //    // for which there does not exist common nodes in the transitions of v and input X

    //    // vertices v in existential vertices
    //    std::vector<Vertex<UBDD>> c1;
    //    for (const auto& v : arena->graph.existential_vertices) {
    //        if (vector_intersection(v.transitions(), X).size() == 0)
    //            c1.push_back(v);
    //    }

    //    // c2(X): vertices v in universal vertices,
    //    // for which the set of transitions from v is a subset of X

    //    // vertices v in universal vertices
    //    std::vector<Vertex<UBDD>> c2;
    //    for (const auto& v : arena->graph.universal_vertices) {
    //        if (subset(v.transitions(), X))
    //            c2.push_back(v);
    //    }

    //    // CPre(X): union between c1 and c2
    //    return vector_union(c1, c2);
    //}


    // Standard algorithm for powerset generation: https://www.geeksforgeeks.org/power-set/
    inline std::vector<std::vector<bool>> powerset(size_t colors) {
        // remake this so that it works with new evaluation function (std::vector<bool>)
        size_t ps_size = pow(2, colors);
        std::vector<std::vector<bool>> ps;
        for (size_t count = 0; count < ps_size; ++count) {
            std::vector<bool> tmp(colors, 0);
            for (size_t bit = 0; bit < colors; ++bit) {
                if (count & (1 << bit))
                    tmp[bit] = true;
            }
            ps.push_back(tmp);
        }
        return ps;
    }

    inline std::vector<size_t> preprocess_to_UBDD(const std::vector<bool>& label) {
        std::vector<size_t> preprocessed;
        for (size_t i = 0; i < label.size(); ++i){
            if (label[i]) preprocessed.push_back(i);
        } return preprocessed;
    }
} // namespace EMHelpers
