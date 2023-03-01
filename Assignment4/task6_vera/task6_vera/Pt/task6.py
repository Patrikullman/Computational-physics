from gpaw import GPAW, restart
from ase import Atoms
from ase.build import fcc111, add_adsorbate
from ase.constraints import FixAtoms
from ase.optimize import BFGS
from ase.visualize import view
from ase.io import read

#reminder: here we should use O, not O2. 

atoms_CO = read('Pt_CO.xyz') #read the files that we got when we did task 3, use ase gui traj, insert CO molecule and saved as xyz file. 
atoms_O = read('Pt_O.xyz')

calc = GPAW(xc='PBE', kpts=(4, 4, 1), parallel={'band': 1}, txt='output_Pt.txt')
atoms_O.set_calculator(calc)
atoms_CO.set_calculator(calc)


relax_Au_CO = BFGS(atoms_CO, trajectory='relax_Pt_CO.traj')
relax_Au_O = BFGS(atoms_O, trajectory='relax_Pt_O.traj')
relax_Au_CO.run(fmax=0.1)  #Relax until the forces acting on each atom are all below 0.1 eV/Å,  fmax=0.01 takes care of this.
relax_Au_O.run(fmax=0.1)

E_surface = -162.4006 #from task 3,  energy of Au  without the adsorbate, this is needed because we need to subtract it from the adsobate energy below.
#but i wonder, do we need the system energy or surface energy from task3..... ?

#energy of the slab with the O atom adsorbed
atoms_O.set_calculator(calc)
E_O_energy = atoms_O.get_potential_energy()  #here we get the energy when we have adsorbed a O atom. 

atoms_CO.set_calculator(calc)
E_CO_energy = atoms_CO.get_potential_energy()

# Calculate the adsorption energies  E_ads = E_adsorbate_on_surface - (E_surface + 1/2 * E_gas)  E_gas we get from task4, i.e 8.749 and 14.210 eV. 
E_ads_O = E_O_energy - (E_surface + 1/2 * 8.749)     
E_ads_CO = E_CO_energy - (E_surface + 1/2 * 14.210)

#now we can use these two E_ads_O and E_ads_CO to get the activation energy (see equation 14 in assignment task).
E_activation = -0.3*(E_ads_O + E_ads_CO) + 0.22 # see equation 14 in A4.
print('Adsorption energy of O atom:', E_ads_O, 'eV')
print('Adsorption energy of CO molecule:', E_ads_CO, 'eV')

print('Activation energy E_activation', E_activation, 'eV')


