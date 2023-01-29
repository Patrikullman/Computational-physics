import gpaw as gp
from ase import Atoms
from ase.optimize import BFGS
from ase.visualize import view
import numpy as np
from gpaw import GPAW, PW

# Create an Atoms object representing the atom structure (this is an ASE object)
atoms = Atoms('Cu2', positions=[[2, 3, 1],[0,1,1]]) #custom-made should maybe not be copper but instead "X", but it doesnt work for some reason. #cell=(10, 10, 10)
atoms.center(vacuum=2.5)
# Set up GPAW calculator
calc = gp.GPAW(mode='fd', h=0.2) #we can play around with different setting on the GPAW calculator to see which one is best.  (fd = finite differences method)
#calc = gp.GPAW(xc='PBE',
#               mode=PW(300),
#               txt='h2.txt')

atoms.set_calculator(calc) #attach calculator to atom, GPAW is a calculator that attaches to an ASE object such as atoms. 

# BFGS optimizes the cluster structure, i.e finds the local minimum (not global minimum)
opt = BFGS(atoms, logfile='bfgs.log')
opt.run(fmax=0.05)

# The cluster is now relaxed and we extract the energy.
energy = atoms.get_potential_energy()
print('Energy: ', energy)

# Save the wavefunction in a .gpw file
calc.write('wavefunction.gpw')

view(atoms)  #plot the cluster. 

