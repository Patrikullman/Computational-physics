#this code comes from the assignment 3 description.
from ase.io import read
from ase import Atoms
from ase.visualize import view
from gpaw import GPAW
from ase.md.npt import NPT
from ase.io import Trajectory 
#a good place to look : https://wiki.fysik.dtu.dk/ase/ase/md.html#   (ase documentation under Molecular dynamics)
atoms = read('test.xyz') #here ill put the snapshot sodium ion in water structure xyz file. 

calc = GPAW(
	#...
	mode = 'lcao',
	xc = 'PBE',
	basis = 'dzp',
	symmetry= {'point_group': False}, # Turn off point -group symmetry
	charge = 1, # Charged system
	txt = 'output.gpaw -out', # Redirects calculator output to this file!
 	)

atoms.set_calculator(calc) #attach calculator to atoms object

from ase.units import fs, kB
dyn = NPT( # Some MD method
	#...
	atoms ,
	pfactor=None,
	temperature_K = 350,
	timestep = 5*fs, # This is not an appropriate timestep , I can tell you that!    (we should use between 1-10 fs) (how can we be sure?)
	ttime = 20*fs, # Don’t forget the fs! (Characteristic timescale of the thermostat)
	externalstress = 0, # We don’t use the barostat , but this needs to be set anyway!
	logfile = 'mdOutput.log', # Outputs temperature (and more) to file at each timestep
	)
trajectory = Trajectory('someDynamics.traj', 'w', atoms)
dyn.attach(trajectory.write , interval =1) # Write the current positions etc. to file each


dyn.run (10) # Run 10 steps of MD simulation
