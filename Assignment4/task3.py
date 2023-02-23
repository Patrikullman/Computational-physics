from ase.build import fcc100
from ase.visualize import view
from ase.build import fcc111
import numpy as np
from ase.build import bulk
from ase import Atoms
from ase.calculators.emt import EMT
from ase.build import fcc111
from gpaw import GPAW, PW
from ase.io import Trajectory
from ase.optimize import BFGS
#https://wiki.fysik.dtu.dk/ase/gettingstarted/surface.html#atoms good link
#https://wiki.fysik.dtu.dk/ase/gettingstarted/tut01_molecule/molecule.html 


#(111)slab of Au,Pt,Rh : use the lattice parameter obtained in task1


# Lattice constants for Au, Pt and Rh (in Å)
a_Au = 4.17 
a_Pt = 3.94
a_Rh = 3.83

# Create (111) surfaces of Au, Pt and Rh
Au = fcc111('Au', a=a_Au, size=(3, 3, 3), vacuum=6.0)
Pt = fcc111('Pt', a=a_Pt, size=(3, 3, 3), vacuum=6.0)
Rh = fcc111('Rh', a=a_Rh, size=(3, 3, 3), vacuum=6.0)

Au.set_pbc((True, True, True))
Pt.set_pbc((True, True, True))  #had to put these to get the pbc, so i can use kpts in GPAW otherwise it complains about pbc. 
Rh.set_pbc((True, True, True))

# Set the distance between periodic slabs in the z-direction to 12 Å
Au.center(vacuum=12.0, axis=2)
Pt.center(vacuum=12.0, axis=2) #leaving 12 Å of empty space around each slab 
Rh.center(vacuum=12.0, axis=2)
#surface_area = Au.get_volume() / Au.get_all_distances().max()
#print(surface_area)
calc1 = GPAW(xc = 'PBE', mode=PW(450),kpts=(4,4,1), txt='calculationAu.txt')   #perform DFT to get energies.  
calc2 = GPAW(xc = 'PBE', mode=PW(450),kpts=(4,4,1), txt='calculationPt.txt')
calc3 = GPAW(xc = 'PBE', mode=PW(450),kpts=(4,4,1), txt='calculationRh.txt')


# Set the calculator for each of the surfaces
Au.set_calculator(calc1)
Pt.set_calculator(calc2)  #this is the lines that are computationally costly, they initiate the DFT calculation by calling the calculator we set above. 
Rh.set_calculator(calc3)

dyn1 = BFGS(Au, trajectory = 'Au_trajectory.traj')
dyn1.run(fmax = 0.05)

dyn2 = BFGS(Pt, trajectory = 'Pt_trajectory.traj')
dyn2.run(fmax = 0.05)

dyn3 = BFGS(Rh, trajectory = 'Rh_trajectory.traj')
dyn3.run(fmax = 0.05)

#trajectory_Au = Trajectory('Au_trajectory.traj', 'w', Au)
#trajectory_Pt = Trajectory('Pt_trajectory.traj', 'w', Pt)
#trajectory_Rh = Trajectory('Rh_trajectory.traj', 'w', Rh)

E_slab_Au = Au.get_potential_energy() #total energy of the whole slab, we get it from the GPAW DFT calculation.
E_slab_Pt = Pt.get_potential_energy() 
E_slab_Rh = Rh.get_potential_energy() 

calc1.write('output_Au.gpw', mode = 'all')
calc2.write('output_Pt.gpw', mode = 'all')
calc3.write('output_Rh.gpw', mode = 'all')

with open('potential_energies.txt', 'w') as fp1:
    fp1.write('Au potential energy: {:.4f} eV\n'.format(E_slab_Au))
    fp1.write('Pt potential energy: {:.4f} eV\n'.format(E_slab_Pt))
    fp1.write('Rh potential energy: {:.4f} eV\n'.format(E_slab_Rh))
fp1.close() 


E_bulk_Au = -3.146
E_bulk_Pt = -6.431
E_bulk_Rh = -7.306

surface_area_Au = Au.get_volume() / Au.get_all_distances().max() #surface area of gold metal (volume/distance gives length ^2 which is area)
surface_area_Pt = Pt.get_volume() / Pt.get_all_distances().max()
surface_area_Rh = Rh.get_volume() / Rh.get_all_distances().max()
N = 27
# N = len(Au_bulk) #Number of atoms in slab
#E_bulk = Au_bulk.get_potential_energy() / len(Au_bulk)   #E_bulk is the energy per atom, this is for Au.


surface_energy_Au = (E_slab_Au - N*E_bulk_Au) / (2*surface_area_Au)
surface_energy_Pt = (E_slab_Pt - N*E_bulk_Pt) / (2*surface_area_Pt)
surface_energy_Rh = (E_slab_Rh - N*E_bulk_Rh) / (2*surface_area_Rh)


#print(surface_energy_Au)
#print(surface_energy_Pt)
#print(surface_energy_Rh)

with open('surface_energies.txt', 'w') as fp2:
    fp2.write('Au surface energy: {:.4f} eV\n'.format(surface_energy_Au))
    fp2.write('Pt surface energy: {:.4f} eV\n'.format(surface_energy_Pt))
    fp2.write('Rh surface energy: {:.4f} eV\n'.format(surface_energy_Rh))
    
fp2.close()

#then we just do this for the two other surfaces, Pt and Rh.
