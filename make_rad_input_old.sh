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

./read_stl solids/Solid_below_crackplane_1.stl
mv Vertices.dat V1.dat
mv Faces.dat F1.dat

./read_stl solids/Solid_above_crackplane_1.stl
mv Vertices.dat V2.dat
mv Faces.dat F2.dat

./get_intersection_mesh $grid $grid > solids/overlap.dat

module load matlab
x=`tail -n 1 chosen.solution.txt | awk '{print $1" "$2" "$3}'`
y=`tail -n 1 chosen.solution.txt | awk '{print $4" "$5" "$6}'`
ans=`matlab -nosplash -nodisplay -r "peaks;quit" | grep -A2 "ans" | tail -n 1 | awk '{print $1" "$2}'`
xmax=`echo $ans | awk '{print $1}'`
ymax=`echo $ans | awk '{print $2}'`
echo "Plane1: Max z along [ "$x" ] = "$xmax >> barriers.dat
echo "Plane1: Max z along [ "$y" ] = "$ymax >> barriers.dat

mv solids/overlap.dat solids/overlap1.dat
mv OverlapPlot.png solids/OverlapPlot1.png
mv Slices.png solids/Slices1.png

./read_stl solids/Solid_below_crackplane_2.stl
mv Vertices.dat V1.dat
mv Faces.dat F1.dat

./read_stl solids/Solid_above_crackplane_2.stl
mv Vertices.dat V2.dat
mv Faces.dat F2.dat

./get_intersection_mesh $grid $grid > solids/overlap.dat

module load matlab
x=`tail -n 1 chosen.solution.txt | awk '{print $1" "$2" "$3}'`
y=`tail -n 1 chosen.solution.txt | awk '{print $4" "$5" "$6}'`
ans=`matlab -nosplash -nodisplay -r "peaks;quit" | grep -A2 "ans" | tail -n 1 | awk '{print $1" "$2}'`
xmax=`echo $ans | awk '{print $1}'`
ymax=`echo $ans | awk '{print $2}'`
echo "Plane2: Max z along [ "$x" ] = "$xmax >> barriers.dat
echo "Plane2: Max z along [ "$y" ] = "$ymax >> barriers.dat

mv solids/overlap.dat solids/overlap2.dat
mv OverlapPlot.png solids/OverlapPlot2.png
mv Slices.png solids/Slices2.png

#cat Assembly_allcracks_1.stl | awk '{ if ($1 == "vertex") {$2=$2+1.0;$3=$3+1.0} print $0 }' > Assembly_allcracks_1.shifted.stl
