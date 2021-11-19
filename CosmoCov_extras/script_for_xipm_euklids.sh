#!/bin/bash
N=8
( 
for thing in {1..903}; do 
   ((i=i%N)); ((i++==0)) && wait
   ./cov_flat_fft $thing ini_files/cov_euklids.ini  &
done
)

