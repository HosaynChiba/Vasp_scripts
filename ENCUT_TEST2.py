import os
import numpy as np
import matplotlib.pyplot as plt
from ase.build import bulk
from ase.spacegroup import crystal
from ase.calculators.vasp import Vasp
a = 3.192
c = 13.378
mos2= crystal(['Mo', 'S'], basis=[(0.3333, 0.3333, 0.75), (0.3333, 0.3333, 0.1331)], spacegroup=194, cellpar=[a, a, c, 90, 90, 120]) 
encuts = [encut for encut in range(200, 840, 50)] # 200~800
energies = []
for index, encut in enumerate(encuts):
    vasp_directory = '{:0>2d}.encut-{}'.format(index, encut)
    if not os.path.exists(vasp_directory):
        os.makedirs(vasp_directory)
    os.chdir(vasp_directory) 
    
    calculator = Vasp(xc='PBE', kpts=(12, 12, 12),
                    encut=encut,
                    ismear=-5,
                    lreal=False, lcharg=False, lwave=False
                    )
    mos2.set_calculator(calculator)
    energy = mos2.get_potential_energy()
    energies.append(energy)
    print '{:20s}: {:10.3f}'.format(vasp_directory, energy)
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
