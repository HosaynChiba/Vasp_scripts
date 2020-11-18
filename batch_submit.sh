#!/bin/bash
path=`pwd`
for i in *.encut-*
do
  pwd
  cd $path/$i
  echo "Getting in to $i"
  #awk '{sub(/ENCUT=/,"ENCUT=$i[-*]")}1' INCAR
  #sed -i "s/ENCUT  =  [0-9][0-9]0/ENCUT  =  ${i:(-3)}/g" INCAR
  #! your vasp job script or vasp_std
  #echo "${i:(-3)}"
  mpirun -np 4 vasp_std
  E=`grep "TOTEN" OUTCAR | tail -1 | awk '{printf "%12.6f \n", $5 }'`
  echo "${i:(-3)}" $E >> ../E-v-ENCUT.dat
  #cp CONTCAR POSCAR
  echo "Exiting $i"
done
shutdown
