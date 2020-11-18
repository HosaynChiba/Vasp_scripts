#!/bin/bash
path=`pwd`
for i in lattice_*
do
  pwd
  cd $path/$i
  echo "Getting in to $i"
  #! your vasp job script or vasp_std
  mpirun -np 4 vasp_std
  #cp CONTCAR POSCAR
  echo "Exiting $i"
done
