#!/usr/bin/env bash

./construct_evaluator.py $1
[ $? -ne 0 ] && exit

make clean
make

time ./bin $2
