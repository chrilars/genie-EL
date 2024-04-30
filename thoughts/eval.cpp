#include "../include/lib/condition_evaluation.hh"

// Test file for condition_evaluation.hh
int main() {
    // Test cases
    std::string test1 = "(!0 | (1)) & ((!0) | (!3)) & (2)";
    std::string test2 = "((!0&1)|(!2&3))&(!4|5)&(!6|7)";
    std::string test3 = "(!0 | 1) & (!2 | 3) & (!4 | 5) & (!6 | 7)";
    std::string test4 = "(!0 | 1 | 2 | 3 | 4 | 5 | 6 | 7) & (!2 | 3 | 4 | 5 | 6 | 7) & (!4 | 5 | 6 | 7) & (!6 | 7)";
    std::string test5 = "(Fin 0 | Inf 1 | Inf 2 | Inf 3 | Inf 4 | Inf 5 | Inf 6 | Inf 7) & (Fin 2 | Inf 3 | Inf 4 | Inf 5 | Inf 6 | Inf 7) & (Fin 4 | Inf 5 | Inf 6 | Inf 7) & (Fin 6 | Inf 7)";
    std::string test6 = "!0 & (1 | 2)";
    std::string test7 = "(!0|1|2|3|4|5|6|7|8)&(!2|3|4|5|6|7|8)&(!4|5|6|7|8)&(!6|7|8)&(!8|9)";
    std::string test8 = "(0&(!1)&(!2)&(!3)&(!4)&(!5)&(!6)&(!7)&(!8))|(2&(!3)&(!4)&(!5)&(!6)&(!7)&(!8))|(4&(!5)&(!6)&(!7)&(!8))|(6&(!7)&(!8))|(8)";
    std::string test9 = "(0 & !1 & !2 & !3 & !4 & !5 & !6 & !7 & !8) | (2 & !3 & !4 & !5 & !6 & !7 & !8) | (4 & !5 & !6 & !7 & !8) | (6 & !7 & !8) | (8)";

    std::vector<std::string> s = tokenize(test9);
    std::vector<std::string> s2 = infix2postfix(s);

    print_tokens(s);
    print_tokens(s2);

    std::vector<bool> colors = {1, 1, 1, 1, 1, 1, 1, 1, 1};
    std::cout << eval_postfix(s2, colors) << std::endl;

    return 0;
}
