#include "condition_evaluation.hh"


// INPUT ALPHABET
// a | !a | (a) | a&a | a|a
// (a = Inf(a), !a = Fin(a))
//
// INPUT EXAMPLE
// 0 & !1 | (1 | 2)
// Numbers represent variable indices, 0 -> variable 0.

void printTokens(std::vector<std::string> input);
std::vector<std::string> tokenize(std::string input);
std::vector<std::string> infix2postfix(std::vector<std::string> tokens);
bool postfixEval (std::vector<std::string> postfix, std::vector<bool> colors);

int main() {
    // Test cases
    std::string test1 = "(!0 | (1)) & ((!0) | (!3)) & (2)";
    std::string test2 = "((!0&1)|(!2&3))&(!4|5)&(!6|7)";
    std::string test3 = "(!0 | 1) & (!2 | 3) & (!4 | 5) & (!6 | 7)";
    std::string test4 = "(!0 | 1 | 2 | 3 | 4 | 5 | 6 | 7) & (!2 | 3 | 4 | 5 | 6 | 7) & (!4 | 5 | 6 | 7) & (!6 | 7)";
    std::string test5 = "(Fin 0 | Inf 1 | Inf 2 | Inf 3 | Inf 4 | Inf 5 | Inf 6 | Inf 7) & (Fin 2 | Inf 3 | Inf 4 | Inf 5 | Inf 6 | Inf 7) & (Fin 4 | Inf 5 | Inf 6 | Inf 7) & (Fin 6 | Inf 7)";
    std::string test6 = "!0 & (1 | 2)";

    std::vector<std::string> s = tokenize(test6);
    std::vector<std::string> s2 = infix2postfix(s);

    printTokens(s);
    //printTokens(s2);
    std::vector<bool> colors = {false, false, true, false};

    std::cout << postfixEval(s2, colors) << std::endl;

    return 0;
}
