#include "ZielonkaTree.hh"
#include "structures.hh"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>
#include <queue>

// Private
void ZielonkaTree::generate() {
    // 1. get colors of root (label)
    // 2. for each subset in powerset(label) which is of opposite winning status as current,
    //    generate a child node with that subset of colors.
    //    This includes evaluating the winning status of the subsets
    // 3. continue while there are nodes which have possible child nodes to be created

    std::queue<ZielonkaNode*> q;
    q.push(root);
    int32_t hierarchy = 1;
    while (!q.empty()) {
        ZielonkaNode* current = q.front();
        q.pop();
        std::vector<std::vector<Color>> ps = powerset(current->label);
        auto compare = [](const std::vector<Color>& v1, const std::vector<Color>& v2) { return v1.size() > v2.size(); };
        std::sort(ps.begin(), ps.end(), compare);

        for (const auto& color_set : ps) {
            // TODO, subsets with smallest change, must be same amount of difference
            // TODO, child nodes from the same parent node cannot be subsets of eachother

            if (eval_condition_dependent(color_set, current->winning)) {
                ZielonkaNode *child_zn = new ZielonkaNode {
                    .parent = current,
                    .label  = color_set,
                    .level  = current->level + 1,
                    .hierarchy = hierarchy++,
                    .winning = !current->winning
                };
                current->children.push_back(child_zn);
                q.push(child_zn);
            }
        }
    }
}


bool ZielonkaTree::eval_condition_dependent(std::vector<Color> colors, bool winning) {
    // TODO, compare colors to phi but with respect to if it is supposed to be a winning node or not
}


bool ZielonkaTree::eval_condition(std::vector<Color> colors) {
    // TODO, just compare colors to phi
}


// Public
ZielonkaTree::ZielonkaTree(Arena *A) {
    // initialize root node from the colors in the arena
    arena = A;
    phi   = A->phi;
    root  = new ZielonkaNode {
        .parent = nullptr,
        .label  = A->graph.colors,
        .level  = 0,
        .hierarchy = 0,
        .winning = eval_condition(A->graph.colors)
    };
}


bool ZielonkaTree::subset(const std::vector<Color>& gamma, const std::vector<Color>& label) {
    // Check if gamma is subset of label
    for (const auto& g : gamma) {
        if (std::find(label.begin(), label.end(), g) == label.end())
            return false;
    }
    return true;
}


std::vector<Vertex> ZielonkaTree::gamma_inverse(std::vector<Color> colors) {
    // return vertices in graph for which gamma is a subset of given set of colors
    std::vector<Vertex> vertices;
    for (const auto& v : arena->graph.vertices) {
        if (subset(v.gamma, colors))
            vertices.push_back(v);
    }
    return vertices;
}


ZielonkaNode ZielonkaTree::anc(ZielonkaNode s, ZielonkaNode t, std::vector<Vertex> vs) {
    // TODO
}

std::vector<Vertex> ZielonkaTree::CPre(std::vector<Vertex> X) {
    // TODO, Controlled predecessor function
    // CPre(X) = {v ∈ V∃ | E(v) ∩ X = ∅} ∪ {v ∈ V∀ | E(v) ⊆ X}
    // c1 = {v ∈ V∃ | E(v) ∩ X = ∅}
    // c2 = {v ∈ V∀ | E(v) ⊆ X}

    // c1(X): vertices v in existential vertices,
    // for which there does not exist common nodes in the transitions of v and input X

    // vertices v in existential vertices
    std::vector<Vertex> c1;
    for (const auto& v : arena->graph.existential_vertices) {
        if (intersection(v.transitions(), X).size() == 0)
            c1.push_back(v);
    }

    // c2(X): vertices v in universal vertices,
    // for which the set of transitions from v is a subset of X

    // vertices v in universal vertices
    std::vector<Vertex> c2;
    for (const auto& v : arena->graph.universal_vertices) {
        if (subset(v.transitions(), X))
            c2.push_back(v);
    }

    // CPre(X): union between c1 and c2
    return union(c1, c2);
}


// Standard algorithm for powerset generation: https://www.geeksforgeeks.org/power-set/
std::vector<std::vector<Color>> ZielonkaTree::powerset(std::vector<Color> colors) {
    size_t size = colors.size();
    size_t ps_size = pow(2, size);
    std::vector<std::vector<Color>> ps(ps_size);
    for (size_t count = 0; count < ps_size; ++count) {
        std::vector<Color> tmp;
        for (size_t bit = 0; bit < size; ++bit) {
            if (count & (1 << bit))
                tmp.push_back(colors[bit]);
        }
        ps.push_back(tmp);
    }
    return ps;
}

