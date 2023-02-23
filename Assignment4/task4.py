from ase.visualize import view
import numpy as np
from ase import Atoms
from gpaw import GPAW, PW
from ase.build import molecule




################This is task 4 AND 5.

#https://wiki.fysik.dtu.dk/ase/tutorials/db/db.html  link on how to calculate adsorption energy. 

'''
CO = molecule('CO', vacuum=12.0) #defining molecule objects, vacuum 12 makes it gas phase like. 

#DFT
calc = GPAW(mode=PW(450),
            xc='PBE',
            kpts=(1, 1, 1),  #(1,1,1) is gamma according to ASE documentation. https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html    #gamma=True,
            txt='task4CO.txt')
CO.set_calculator(calc)

E_CO = CO.get_potential_energy()

print("Energy of CO molecule in gas phase:", E_CO, "eV")
view(CO)
'''

O2 = molecule('O2', vacuum=12.0) #defining molecule objects, vacuum 12 makes it gas phase like. 

calc = GPAW(mode=PW(450),
            xc='PBE',
            kpts=(1, 1, 1),  #(1,1,1) is gamma according to ASE documentation. https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html    #gamma=True,
            txt='task4O2.txt')
O2.set_calculator(calc)




################################Task 5#######################


from ase import Atoms
#from ase.calculators.gpaw import GPAW
from ase.vibrations import Vibrations
from gpaw import GPAW, PW
from ase.build import molecule
import numpy as np
from ase.optimize import BFGS
from ase.thermochemistry import IdealGasThermo
from ase.optimize import QuasiNewton

#basically we get the vibrational energies by using ase.vibrations library. We get potential energy by DFT GPAW and we get entropy by using ase.thermchemistry library.

#good links: https://wiki.fysik.dtu.dk/ase/ase/thermochemistry/thermochemistry.html and https://wiki.fysik.dtu.dk/gpaw/tutorialsexercises/vibrational/vibrations/vibrations.html 

potentialenergy = O2.get_potential_energy()
print("Energy of O2 molecule in gas phase:", potentialenergy, "eV")
dyn = QuasiNewton(O2)  #eller lägg innan get_potential energy ? 
dyn.run(fmax=0.05)
vib = Vibrations(O2)
vib.run()
vib.summary()
for mode in range(6):
    vib.write_mode(mode)   #try ase gui vib.5.traj to look at the vibrations. 

vib_energies = vib.get_energies() #collecting vibrational energies.

thermo = IdealGasThermo(vib_energies=vib_energies,          #https://wiki.fysik.dtu.dk/ase/ase/thermochemistry/thermochemistry.html    here you can read about the ase.thermochemistry library.
                        potentialenergy=potentialenergy,
                        geometry='linear',atoms=O2,
                        symmetrynumber=2, spin=1) #should this be 0 or 1 ?

#G = thermo.get_gibbs_energy(temperature=298.15, pressure=101325.) #probably dont need this but i was curious.
E = thermo.get_entropy(temperature=298.15, pressure=101325.)  # entropy extracted by from ase.thermochemistry library.








