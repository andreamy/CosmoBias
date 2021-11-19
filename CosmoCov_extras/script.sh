#!/bin/bash
N=12

for thing in {1..465}; do 
   ((i=i%N)); ((i++==0)) && wait
   ./cov_flat_fft $thing ini_files/ng_cov_kids.ini  &
done


