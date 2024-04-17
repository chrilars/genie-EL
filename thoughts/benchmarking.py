#!/usr/bin/env python3
# Produce benchmarking conditions

import sys


def rabin_condition(n: int) -> str:
    # produce rabin condition of size n
    condition = ""
    for i in range(0, 2*n, 2):
        if i:
            condition += "|"
        condition += f"(!{i} & {i+1})"
    return condition


def other_condition(n: int) -> str:
    # produce the other generalized condition
    condition = ""
    for k in range(0, 2*n, 2):
        if k:
            condition += "&"
        condition += f"(!{k}|"
        for i in range(2*n-k-1):
            if i:
                condition += "|"
            condition += f"{k+i+1}"
        condition += ")"
    return condition


if __name__ == '__main__':
    if len(sys.argv) < 2:
        exit(1)
    rc = rabin_condition(int(sys.argv[1]))
    with open("tmp_rabin.cfg", "w") as rf:
        rf.write(rc)
    oc = other_condition(int(sys.argv[1]))
    with open("tmp_other.cfg", "w") as of:
        of.write(oc)
