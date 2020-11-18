#!/bin/bash
a=0
E=0
V=0
for d in ./relax_*
do  
	cd "$d" 
	echo "Getting in to $d"
	cp INCAR.relax INCAR
	mpirun -np 4 vasp_std
	cp CONTCAR POSCAR
	cp INCAR.static INCAR
	mpirun -np 4 vasp_std
	E=`grep "TOTEN" OUTCAR | tail -1 | awk '{printf "%12.6f \n", $5 }'`
	V=`grep "volume" OUTCAR | tail -1 | awk '{printf "%12.4f \n" , $5}'`
	a=`sed -n '3p' CONTCAR | awk '{printf "%12.4f \n" , $1}'`
	echo $a $V $E >> ../SUMMARY.dat
	#sed -n '1p' SUMMARY.dat >> ../SUMMARY.dat
	echo "Exiting $d"
	cd ..
done
#echo $a $V $E >> SUMMARY.dat
#shutdown
