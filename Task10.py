import gpaw as gp
from ase.io import read
from ase import Atoms
from ase.optimize import BFGS
from ase.visualize import view
import numpy as np
from ase.io import write
from gpaw import GPAW, PW
#,[1, 1, 1],[0,0,1],[1, 0, 1],[0,0,2]
# Create an Atoms object representing the atom structure (this is an ASE object)
atoms = read('lowest_energy_structure_Na6.xyz') # lowest_energy_structure_Na6.xyz ,lowest_energy_structure_Na7.xyz, lowest_energy_structure_Na8.xyz   (here we load the xyz files we got from 6,7,8)


atoms.center(vacuum=3.0) 
# Set up GPAW calculator
#calc = gp.GPAW(mode='fd',xc='LDA', h=0.2,txt='h2.txt') #we can play around with different setting on the GPAW calculator to see which one is best.  (fd = finite differences method)
calc = gp.GPAW(xc='PBE', 
               mode=PW(490),
               txt='h2.txt')

atoms.set_calculator(calc) #attach calculator to atom, GPAW is a calculator that attaches to an ASE object such as atoms. 

# BFGS optimizes the cluster structure, i.e finds the local minimum (not global minimum)
opt = BFGS(atoms, logfile='bfgs.log')
opt.run(fmax=0.05)

# The cluster is now relaxed and we extract the energy.
energy = atoms.get_potential_energy()
print('Energy: ', energy)


# Save the wavefunction in a .gpw file
calc.write('wavefunctions.gpw')

view(atoms)  #plot the cluster. 

basename = 'wavefunction'
#atoms, calc = restart('wavefunction.gpw')

# loop over all wfs and write their cube files  ( see this link : https://wiki.fysik.dtu.dk/gpaw/tutorialsexercises/wavefunctions/plotting/plot_wave_functions.html ) 
nbands = calc.get_number_of_bands() 
for band in range(nbands):
    wf = calc.get_pseudo_wave_function(band=band)
    fname = f'{basename}_{band}.cube'
    print('writing wf', band, 'to file', fname)
    write(fname, atoms, data=wf)  
