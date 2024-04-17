#!/usr/bin/env python3

# INPUT ALPHABET
# a | !a | (a) | a&a | a|a
# (a = Inf(a), !a = Fin(a))
#
# INPUT EXAMPLE
# 0 & !1 | (1 | 2)
# Numbers represent variable indices, 0 -> variable 0.

import re
import sys

OUT_FILE = "condition_evaluator.hh"


def convert(in_file: str):
    with open(in_file, 'r') as ifile:
        formula = ifile.read().strip()
        print("INPUT:")
        print(formula)
        print()
        formula = formula.replace('&', "&&")
        formula = formula.replace('|', "||")
        formula = formula.replace("Fin ", "!")
        formula = formula.replace("Fin", "!")
        formula = formula.replace("Inf ", "")
        formula = formula.replace("Inf", "")

        rgx = r"\d+"
        nums = set(int(n) for n in re.findall(rgx, formula))
        max_nums = max(nums)
        if max_nums >= len(nums):
            raise SyntaxError(f"Variables are not in order!\nMax variable: {max_nums}, number of variables: {len(nums)}.")
        i = 0
        final = ""
        len_formula = len(formula)
        while i < len_formula:
            c = formula[i]
            if not c.isdigit():
                final += c
            else:
                done = False
                final += "vars["
                while i < len_formula and formula[i].isdigit():
                    final += formula[i]
                    if i == len_formula - 1:
                        final += "]"
                        done = True
                    i += 1
                if not done:
                    final += "]"
                continue
            i += 1
        print("OUTPUT:")
        print(final)
        print()

    with open(OUT_FILE, 'w') as ofile:
        function_template = "#pragma once\n\n"
        function_template += "#include <vector>\n"
        function_template += "inline bool condition_evaluator(const std::vector<bool>& vars) {\n"
        function_template += f"    return {final};\n"
        function_template += "}\n"
        ofile.write(function_template)
        print("FUNCTION_TEMPLATE:")
        print(function_template)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        raise RuntimeError("No formula file provided!")
    in_file = sys.argv[1]
    convert(in_file)
