#!/bin/bash
N=8

( 
for thing in {4 5 7 8 10 27 41 45 48 49 50 51 53 54}; do 
   ((i=i%N)); ((i++==0)) && wait
   ./cov_flat_fft $thing chi2tests/test_$1.ini  &
done
)

