import os
import numpy as np
import matplotlib.pyplot as plt
from ase.build import bulk
from ase.io import read, write
from ase.atoms import Atoms
from ase.calculators.vasp import Vasp
from ase.build import mx2
#mos2 = read('POSCAR')
#metal_symbol = 'Pt'
#structure = 'fcc'
#a0 = 3.98
#metal = bulk(metal_symbol, structure, a=a0)
atoms_symbol = "MoS2"
atoms_kind = "2H"
a_guess = 12.76
atoms = mx2(atoms_symbol, atoms_kind, a=a_guess, thickness=3.19, vacuum=15)
encuts = [encut for encut in range(200, 840, 50)] # 200~800
energies = []
for index, encut in enumerate(encuts):
    vasp_directory = '{:0>2d}.encut-{}'.format(index, encut)
    if not os.path.exists(vasp_directory):
        os.makedirs(vasp_directory)
    os.chdir(vasp_directory) 
    calculator = Vasp(xc='PBE', kpts=(12, 12, 1),
                    encut=encut,
                    ismear=-5,
                    lreal=False, lcharg=False, lwave=False
                    )
    atoms.set_calculator(calculator)
    energy = atoms.get_potential_energy()
    energies.append(energy)
    print('{:20s}: {:10.3f}'.format(vasp_directory, energy))
    os.chdir(os.path.dirname(os.getcwd())) # back to root directory
plt.plot(encuts, energies, '-o')
plt.xlabel('Cutoff Energy (eV)')
plt.ylabel('Total Energy (eV)')
plt.title('PBE encut test')
plt.savefig('encut.png')
plt.close()   
'''
00.encut-200          :     -6.024
01.encut-250          :     -6.103
02.encut-300          :     -6.101
03.encut-350          :     -6.096
04.encut-400          :     -6.096
05.encut-450          :     -6.096
06.encut-500          :     -6.097
07.encut-550          :     -6.097
08.encut-600          :     -6.097
09.encut-650          :     -6.097
10.encut-700          :     -6.097
11.encut-750          :     -6.097
12.encut-800          :     -6.097
'''
