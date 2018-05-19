#!/bin/sh
clear
echo 'The Tomo Program is start!'
echo 'Please ensure the parfile is done:'

gcc -o step1 step_self_new_1.c -lm -w
./step1

gcc -o step2 step_self_new_2.c -lm -w
./step2

gcc -o step3 step_self_new_3.c -lm -w
./step3

gcc -o step4 step_self_new_4.c -lm -w
./step4

rm *~
