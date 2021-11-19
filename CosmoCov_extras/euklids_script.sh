#!/bin/bash
N=12
 
for thing in {1..465}; do 
   ((i=i%N)); ((i++==0)) && wait
   ./cov $thing ini_files/cov_euklids.ini  &
done


