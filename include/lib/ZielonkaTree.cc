#include "ZielonkaTree.hh"
#include "ELHelpers.hh"
#include "structures.hh"
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <queue>
#include "condition_evaluator.hh"


bool cmp_descending_count_true(const std::vector<bool>& a, const std::vector<bool>& b) {
    return std::count(a.begin(), a.end(), true) > std::count(b.begin(), b.end(), true);
}

// Private
size_t leaves = 0;
size_t nodes = 0;
void ZielonkaTree::generate() {
    //std::cout << "generating... \n";
    std::queue<ZielonkaNode*> q;
    q.push(root);
    std::vector<std::vector<bool>> ps = ELHelpers::powerset(root->label.size());
    std::sort(ps.begin(), ps.end(), cmp_descending_count_true);
    size_t order = root->order + 1;
    std::vector<std::vector<bool>> seen_from_parent{};
    while (!q.empty()) {
        seen_from_parent.clear();
        ZielonkaNode* current = q.front();
        q.pop();
        nodes++;
        for (size_t i = 0; i < ps.size(); ++i) {
            std::vector<bool> color_set = ps[i];
            if (!ELHelpers::proper_subset(color_set, current->label))
                continue;
            bool seen = false;
            for (const auto& s : seen_from_parent) {
                if (ELHelpers::proper_subset(color_set, s)) {
                    seen = true;
                    break;
                }
            }
            if (seen) continue;
            if (eval_condition(color_set) != current->winning) {
                ZielonkaNode *child_zn = new ZielonkaNode {
                    .children = {},
                    .child_differences = {},
                    .parent = current,
                    .parent_order = current->order,
                    .label = color_set,
                    .level = current->level + 1,
                    .order = order++,
                    .winning = !(current->winning)
                };
                current->child_differences.push_back(
                    ELHelpers::label_difference(current->label, color_set)
                );
                seen_from_parent.push_back(color_set);
                current->children.push_back(child_zn);
                q.push(child_zn);
            }
        }
        if (current->children.empty()) leaves++;
    }
    //std::cout << "done generating...\n";
}
void ZielonkaTree::generate_parity() {
    // firstly evaluate root, then remove the last color from the current colorset
    //std::cout << "generating... \n";
    size_t order = root->order + 1;
    std::vector<bool> colors;
    for (const bool& b : root->label)
        colors.push_back(b);
    ZielonkaNode* current = root;
    for (int i = colors.size()-1; i >= 0; --i) {
        colors[i] = false;
        ZielonkaNode *child_zn = new ZielonkaNode {
            .children = {},
            .child_differences = {},
            .parent = current,
            .parent_order = current->order,
            .label = colors,
            .level = current->level + 1,
            .order = order++,
            .winning = !(current->winning)
        };
        current->child_differences.push_back(
            ELHelpers::label_difference(current->label, colors)
        );
        current->children.push_back(child_zn);
        current = child_zn;
    }
}

std::string label_to_string(std::vector<bool> label) {
    std::string s;
    for (size_t i = 0; i < label.size(); ++i) {
        if (label[i])
            s += static_cast<char>('a' + i);
    }
    if (s.empty()) return "∅";
    return s;
}

bool ZielonkaTree::eval_condition(std::vector<bool> colors) {
    return condition_evaluator(colors);
}

void ZielonkaTree::displayZielonkaTree() {
    // BFS
    std::cout << "displaying...\n";
    std::queue<ZielonkaNode*> q;
    q.push(root);
    while (!q.empty()) {
        ZielonkaNode* current = q.front();
        ELHelpers::print_label(current);
        std::cout << "from: "
                  << current->parent_order
                  << ", order: "
                  << current->order << '\n';
        for (size_t i = 0; i < current->children.size(); ++i) {
            std::cout << label_to_string(current->child_differences[i]) << '\n';
        }
        std::cout << '\n';
        q.pop();
        for (ZielonkaNode *z : current->children) {
            q.push(z);
        }
    }
}




void printNTree(ZielonkaNode* x, std::vector<bool> flag, int depth = 0, bool isLast = false) {
    // Taken from https://www.geeksforgeeks.org/print-n-ary-tree-graphically/
    // Condition when node is None
    if (x == NULL)
        return;

    // Loop to print the depths of the
    // current node
    for (int i = 1; i < depth; ++i) {
        // Condition when the depth
        // is exploring
        if (flag[i] == true) {
            std::cout << "│ "
                << " "
                << " "
                << " ";
        }
        // Otherwise print 
        // the blank spaces
        else {
            std::cout << " "
                << " "
                << " "
                << " ";
        }
    }
    // Condition when the current
    // node is the root node
    if (depth == 0)
        std::cout << label_to_string(x->label) << " " << (x->winning? 'W' : 'L') << '\n';

    // Condition when the node is 
    // the last node of 
    // the exploring depth
    else if (isLast) {
        std::cout << "└── " << label_to_string(x->label) << " " << (x->winning? 'W' : 'L') << '\n';
         
        // No more childrens turn it 
        // to the non-exploring depth
        flag[depth] = false;
    }
    else {
        std::cout << "├── " << label_to_string(x->label) << " " << (x->winning? 'W' : 'L') << '\n';
    }
 
    size_t it = 0;
    for (auto i = x->children.begin();
    i != x->children.end(); ++i, ++it)
 
        // Recursive call for the
        // children nodes
        printNTree(*i, flag, depth + 1, 
            it == (x->children.size()) - 1);
    flag[depth] = true;
}

void ZielonkaTree::graphZielonkaTree() {
    printNTree(root, std::vector<bool>(85, true));
}

// Public
ZielonkaTree::ZielonkaTree(size_t colors) {
    std::vector<bool> label(colors, true);
    root = new ZielonkaNode {
        .children  = {},
        .child_differences = {},
        .parent = nullptr,
        .parent_order = 0,
        .label = label,
        .level = 1,
        .order = 1,
        .winning = eval_condition(label)
    };
    //generate();
    generate_parity();
    //graphZielonkaTree();
    std::cout << "leaves: "<< leaves << '\n';
    std::cout << "nodes: " << nodes  << '\n';
    //displayZielonkaTree();
    //graphZielonkaTree();
}




//ZielonkaTree::ZielonkaTree(Arena *A) {
//    // initialize root node from the colors in the arena
//    arena = A;
//    phi   = A->phi;
//    root  = new ZielonkaNode {
//        .parent = nullptr,
//        .label  = A->graph.colors,
//        .level  = 0,
//        .order = 0,
//        .winning = eval_condition(A->graph.colors)
//    };
//}


ZielonkaNode* ZielonkaTree::get_root() { return root; }


ZielonkaNode* ZielonkaTree::leads_to(ZielonkaNode *s, ZielonkaNode *t) {
    // Find s_t: child of s which leads to t
    // (bottom up) start at t and go through parents until s is found
    // then take the last visited
    // TODO, Only works if they are on a straight path, not sure if this is the case (might not be actually)
    size_t low_level = s->level;
    size_t target_order = s->order;
    ZielonkaNode *last_visited = t;
    ZielonkaNode *current = t;
    while (current->level >= low_level) {
        if (current->order == target_order) {
            return last_visited;
        } else {
            last_visited = current;
            if (current->parent == nullptr) {
                throw std::runtime_error("Could not find parent ZielonkaNode, gone above root!");
            }
            current = current->parent;
        }
    }
    throw std::runtime_error("Could not find parent ZielonkaNode, lowest level surpassed!");
}


//template <class UBDD>
//UBDD ZielonkaTree::anc(ZielonkaNode *s, ZielonkaNode *t) {
//    // v1: gamma inverse, subset of l(s)
//    // v2: gamma inverse, not subset of l(s_t)
//    // return intersection(v1, v2)
//    if (s->order == t->order)
//        return ELHelpers::gamma_inverse(t->label, arena, true);
//    ZielonkaNode *s_t = leads_to(s, t);
//    UBDD v1 = ELHelpers::gamma_inverse(  s->label, arena, true);
//    UBDD v2 = ELHelpers::gamma_inverse(s_t->label, arena, false);
//    return ELHelpers::vector_intersection(v1, v2);
//}
//
//
//std::vector<Vertex> ZielonkaTree::solve(ZielonkaNode *s, std::vector<std::vector<Vertex>> vs) {
//    std::vector<std::vector<Vertex>> Xs;
//    std::vector<std::vector<Vertex>> W;
//    std::vector<Vertex> U;
//    if (!s->winning) Xs = vs;                    // if s ∈ T_circle then Xs ← ∅ else Xs ← V
//    W = ELHelpers::vector_difference(vs, Xs);    // W ← V \ Xs
//    while (Xs != W) {                                               // while Xs != W do
//        Xs = W;                                                     // Xs ← W
//        if (!s->children.empty()) {                                 // if R(s) != ∅ then
//            for (ZielonkaNode *t : s->children) {                  // for t ∈ R(s) do
//                U = solve(t, W);
//                if (s->winning)
//                    Xs = ELHelpers::vector_union(Xs, U);
//                else
//                    Xs = ELHelpers::vector_intersection(Xs, U);
//            }
//        } else {
//            // TODO, do whatever the solver wants
//        }
//    }
//}

