#pragma once

#include <cstddef>
#include <vector>

struct ZielonkaNode {
    std::vector<ZielonkaNode*> children;
    std::vector<std::vector<bool>> child_differences;
    ZielonkaNode *parent;
    size_t parent_order;
    std::vector<bool> label;
    size_t level;
    size_t order;
    bool winning;
};
