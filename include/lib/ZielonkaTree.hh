#pragma once

#include "structures.hh"
#include <vector>
#include "condition_evaluation.hh"


class ZielonkaTree {
private:
    // Private Variables
    ZielonkaNode *root;
    std::vector<std::string> phi; // Emerson-Lei condition in tokenized postfix format

    // Private methods
    void generate();
    void generate_parity();
    void generate_phi(std::string);
    bool eval_condition(std::vector<bool>);
    void displayZielonkaTree();
    void graphZielonkaTree();

public:
    ZielonkaTree(std::string, size_t);
    ~ZielonkaTree() {};

    ZielonkaNode* get_root();
    ZielonkaNode* leads_to(ZielonkaNode*, ZielonkaNode*);
};

