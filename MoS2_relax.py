import numpy as np
import matplotlib
matplotlib.use('Agg')
import os
from ase.build import mx2
from ase.calculators.vasp import Vasp
from ase.eos import EquationOfState
atoms_symbol = "MoS2"
atoms_kind = "2H"
a_guess = 3.12
atoms = mx2(atoms_symbol, atoms_kind, a=a_guess, thickness=3.19, vacuum=7.5)
h = atoms.get_cell_lengths_and_angles()[2]
print 'MoS2-cell-c : {:.3f}'.format(h)
lattice_a = np.linspace(a_guess*0.95, a_guess*1.05, 11)
volumes = []
energies = []
for index, a in enumerate(lattice_a):
    vasp_directory = '{:0>2d}.lattice-{:.3f}'.format(index, a)
    if not os.path.exists(vasp_directory):
        os.makedirs(vasp_directory)
    os.chdir(vasp_directory) 
    atoms = mx2(atoms_symbol, atoms_kind, a=a, thickness=3.19, vacuum=7.5)
    calculator = Vasp(xc='PBE', kpts=(6, 6, 1), encut=400,      
                    ismear=0, sigma=0.05,
                    ibrion=1, nsw=200,
                    lreal=False, lcharg=False, lwave=False
                    )
    atoms.set_calculator(calculator)
    volume = atoms.get_volume()
    energy = atoms.get_potential_energy()
    volumes.append(volume)
    energies.append(energy)
    print '{:20s}: {:10.3f} {:10.3f}'.format(vasp_directory, volume, energy)
    os.chdir(os.path.dirname(os.getcwd())) # back to root directory
    
eos = EquationOfState(volumes, energies)
v0, e0, B = eos.fit()
print "v0 = {:.3f}".format(v0)
a0 = np.sqrt( (2*v0) / (h*np.sqrt(3)))
print "a0 = {:.3f}".format(a0)
eos.plot('{}-{}-eos.png'.format(atoms_symbol,atoms_kind))
print "Relaxtion after getting a0"
vasp_directory = 'opt'
if not os.path.exists(vasp_directory):
    os.makedirs(vasp_directory)
os.chdir(vasp_directory) 
atoms = mx2(atoms_symbol, atoms_kind, a=a0, thickness=3.19, vacuum=7.5)
calculator = Vasp(xc='PBE', kpts=(6, 6, 1), encut=400,
                ismear=0, sigma=0.05,
                ibrion=1, nsw=200,
                lreal=False, lcharg=False, lwave=False
                )
atoms.set_calculator(calculator)
final_energy = atoms.get_potential_energy()
thickness = atoms[1].z - atoms[2].z
print "thickness = {:.3f}".format(thickness)
os.chdir(os.path.dirname(os.getcwd())) # back to root directory
'''
MoS2-cell-c : 18.190
00.lattice-2.964    :    138.395    -21.303
01.lattice-2.995    :    141.324    -21.440
02.lattice-3.026    :    144.283    -21.553
03.lattice-3.058    :    147.274    -21.643
04.lattice-3.089    :    150.294    -21.711
05.lattice-3.120    :    153.346    -21.758
06.lattice-3.151    :    156.428    -21.786
07.lattice-3.182    :    159.541    -21.795
08.lattice-3.214    :    162.685    -21.787
09.lattice-3.245    :    165.859    -21.763
10.lattice-3.276    :    169.064    -21.723
v0 = 159.618
a0 = 3.183
Relaxtion after getting a0
thickness = 3.127
'''

