#pragma once

#include "structures.hh"
#include <vector>

class ZielonkaTree {
private:
    // Private Variables
    ZielonkaNode *root;

    // Private methods
    void generate();
    void generate_parity();
    bool eval_condition(std::vector<bool>);
    void displayZielonkaTree();
    void graphZielonkaTree();

public:
    ZielonkaTree(size_t);
    ~ZielonkaTree() {};

    ZielonkaNode* get_root();
    ZielonkaNode* leads_to(ZielonkaNode*, ZielonkaNode*);
};

