#!/usr/bin/env bash

./benchmarking.py $1

#echo "size $1"
#echo "rabin condition"
#cat tmp_rabin.cfg
#echo
#echo "other condition"
#cat tmp_other.cfg
#echo

#./build.sh tmp_rabin.cfg $(($1 * 2))
./build.sh tmp_other.cfg $(($1 * 2))
