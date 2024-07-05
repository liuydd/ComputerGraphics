#!/usr/bin/env bash

cmake -B build
cmake --build build

mkdir -p output
# build/PA1 testcases/basic.txt output/basic_no_nee_10.bmp 10 pt
# build/PA1 testcases/basic.txt output/basic_nee_10.bmp 10 pt
# build/PA1 testcases/basic.txt output/basic_100_r_nee.bmp 100 pt
# build/PA1 testcases/basic.txt output/basic_1000.bmp 1000 pt

# build/PA1 testcases/dof.txt output/dof_10.bmp 10 pt
# build/PA1 testcases/dof.txt output/dof_100.bmp 100 pt
# build/PA1 testcases/dof_compare.txt output/dof_compare_10.bmp 10 pt
# build/PA1 testcases/dof_compare.txt output/dof_compare_100.bmp 100 pt
# build/PA1 testcases/dof.txt output/dof_1k.bmp 1000 pt
# build/PA1 testcases/dof_compare.txt output/dof_compare_1k.bmp 1000 pt
# build/PA1 testcases/move.txt output/move_10.bmp 10 pt
# build/PA1 testcases/move.txt output/move_100.bmp 100 pt
# build/PA1 testcases/move.txt output/move_1000.bmp 1000 pt

# build/PA1 testcases/aotu.txt output/aotu_bricks_100.bmp 100 pt
# build/PA1 testcases/aotu.txt output/aotu_bricks_1000.bmp 1000 pt

build/PA1 testcases/test.txt output/test_no.bmp 100 pt

# build/PA1 testcases/scene09_norm.txt output/norm_10.bmp 10 pt
# build/PA1 testcases/scene10_wineglass.txt output/wineglass_10.bmp 10 pt

# build/PA1 testcases/vase.txt output/vase_1k.bmp 1000 pt
# build/PA1 testcases/vase.txt output/vase_200.bmp 200 pt
# build/PA1 testcases/vase2.txt output/vase2_500.bmp 500 pt
# build/PA1 testcases/vase3.txt output/vase3_shuiwen.bmp 10 pt
# build/PA1 testcases/aotu.txt output/aotu_glossy_2.bmp 1000 pt