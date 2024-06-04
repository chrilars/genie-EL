#!/usr/bin/env bash

python3 benchmarking.py $1
[ $? -ne 0 ] && exit


./a.out $1
