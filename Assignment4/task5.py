from ase import Atoms
#from ase.calculators.gpaw import GPAW
from ase.vibrations import Vibrations
from gpaw import GPAW, PW
from ase.build import molecule
import numpy as np
from ase.optimize import BFGS
from ase.thermochemistry import IdealGasThermo

#basically we get the vibrational energies by using ase.vibrations library. We get potential energy by DFT GPAW and we get entropy by using ase.thermchemistry library.

#good links: https://wiki.fysik.dtu.dk/ase/ase/thermochemistry/thermochemistry.html and https://wiki.fysik.dtu.dk/gpaw/tutorialsexercises/vibrational/vibrations/vibrations.html 

O2 = molecule('O2', vacuum=12.0)

calc_O2 = GPAW(mode='lcao',symmetry={'point_group': False},basis='dzp', xc='LDA', txt='O2.txt')
O2.set_calculator(calc_O2)
potentialenergy = O2.get_potential_energy()

vib = Vibrations(O2)
vib.run()
vib.summary()
for mode in range(6):
    vib.write_mode(mode)   #try ase gui vib.5.traj to look at the vibrations. 

vib_energies = vib.get_energies() #collecting vibrational energies.

thermo = IdealGasThermo(vib_energies=vib_energies,          #https://wiki.fysik.dtu.dk/ase/ase/thermochemistry/thermochemistry.html    here you can read about the ase.thermochemistry library.
                        potentialenergy=potentialenergy,
                        geometry='linear',atoms=O2,
                        symmetrynumber=2, spin=1)

G = thermo.get_gibbs_energy(temperature=298.15, pressure=101325.) #probably dont need this but i was curious.
E = thermo.get_entropy(temperature=298.15, pressure=101325.)  # entropy extracted by from ase.thermochemistry library.





















'''
# Define the CO molecule
CO = Atoms('CO', positions=[(0, 0, 0), (1.128, 0, 0)], cell=[(12, 0, 0), (0, 12, 0), (0, 0, 12)])
CO.set_pbc((True, True, True))
# Set up the GPAW calculator for CO
calc_CO = GPAW(mode='lcao', xc='PBE', kpts=(8, 8, 8), txt='CO.txt')
CO.set_calculator(calc_CO)

# Calculate the vibrational modes and frequencies for CO
vib_CO = Vibrations(CO, name='CO')
vib_CO.run()

# Compute the vibrational entropy of O2 and CO at 300 K and 1 bar
T = 300 # temperature in K
P = 1e5 # pressure in Pa

S_O2 = vib_O2.get_entropy(T=T, P=P) / 4.184 # convert from J/mol/K to cal/mol/K
S_CO = vib_CO.get_entropy(T=T, P=P) / 4.184 # convert from J/mol/K to cal/mol/K

print("Entropy of O2 at 300 K and 1 bar:", S_O2, "cal/mol/K")
print("Entropy of CO at 300 K and 1 bar:", S_CO, "cal/mol/K")
'''
