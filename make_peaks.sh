#!/bin/bash
module load matlab
x=`tail -n 1 chosen.solution.txt | awk '{print $1" "$2" "$3}'`
y=`tail -n 1 chosen.solution.txt | awk '{print $4" "$5" "$6}'`
xy=`tail -n 1 chosen.solution.txt | awk '{print ($1+$4)" "($2+$5)" "($3+$6)}'`
matlab -nosplash -nodisplay -r "peaks1;quit"
#ans=`matlab -nosplash -nodisplay -r "peaks1;quit" | grep -A2 "ans" | tail -n 1 | awk '{print $1" "$2}'`
