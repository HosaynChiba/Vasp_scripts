from ase.utils.eos import EquationOfState
volumes=[278.3900, 217.0600, 183.1600, 61.4400, 118.0700, 45.7400, 74.4100]
energies=[-43.694690, -43.694591, -43.695367, -22.701505, -43.692637, -7.340815, -39.060444]
eos = EquationOfState(volumes, energies)
#v0, e0, B = eos.fit()
#print(('''
#v0 = {0} A^3
#E0 = {1} eV
#B  = {2} eV/A^3'''.format(v0, e0, B)))

eos.plot('mos2-eos.png')
