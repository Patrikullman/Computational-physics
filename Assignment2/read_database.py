import ase
from ase.db import connect
db = connect('gadb.db')
atoms = db.get('id=31').toatoms()  #id 31 is the id corresponding to the lowest structure. 


ase.io.write('lowest_energy_structure_Na8.xyz', atoms)
print(atoms)
