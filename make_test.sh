#!/bin/bash
module load matlab
x=`tail -n 1 chosen.solution.txt | awk '{print $1" "$2" "$3}'`
y=`tail -n 1 chosen.solution.txt | awk '{print $4" "$5" "$6}'`
ans=`matlab -nosplash -nodisplay -nojvm -r "peaks;quit" | grep -A2 "ans" | tail -n 1 | awk '{print $1" "$2}'`
xmax=`echo $ans | awk '{print $1}'`
ymax=`echo $ans | awk '{print $2}'`
echo "Max z along [ "$x" ] = "$xmax > barriers.dat
echo "Max z along [ "$y" ] = "$ymax >> barriers.dat

