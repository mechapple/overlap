#!/bin/bash
rm out.lmp solids/* barriers.dat
make clean
make
./gccm_generate
rfac=0.1
mol_size=28
grid=10
read nx ny nz <<<$(echo 5 5 3)
atomsk newcell.lmp -duplicate $nx $ny $nz out.lmp
cat out.lmp | awk 'NF==5' | awk -v factor=$rfac 'BEGIN {arr[1]=1.7; arr[2]=1.55; arr[3]=1.52; arr[4]=1.10;} {printf "%d %4.8f %4.8f %4.8f %2.3f\n", $1,$3,$4,$5,factor*arr[$2]}' > hmx.dat
awk -v apm=$mol_size '{sumx+=$2;sumy+=$3;sumz+=$4} (NR%apm)==0 {printf("%d %4.10g %4.10g %4.10g\n",NR/apm,sumx/apm,sumy/apm,sumz/apm); sumx=0;sumy=0;sumz=0;}' hmx.dat > hmx_coms.dat 
grep -A3 "xlo" out.lmp > bounds.dat
grep -A3 "xlo" newcell.lmp > bounds.cell
./radical_hmx
#povray hmx.pov
mkdir -p solids
./create_assembly_stl $mol_size $nx $ny $nz > solids/out
