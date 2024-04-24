#pragma once

#include <vector>
inline bool condition_evaluator(const std::vector<bool>& vars) {
    return (!vars[0] || (vars[1])) && ((!vars[0]) || (!vars[3])) && (vars[2]);
}
