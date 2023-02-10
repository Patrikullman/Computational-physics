import numpy as np
from matplotlib import pyplot as plt
from ase.io.trajectory import Trajectory
from ase import Atoms

atoms = Trajectory('NaCluster24.traj') #someDynamics NaCluster24

distances = [] #np.array([])
for i in range(1000,14000):   #i think we need to remove the first time frames since we need to pick frames after equilibration.
    distances = np.append(distances,atoms[i].get_distances(72,atoms[i].get_atomic_numbers()==8,mic=True))

histogram = np.histogram(distances,bins=1000)

V = atoms[0].get_volume()  #cell volume, this is an ase function. 
N = 24 #oxygen molecules
rho = N/V #number density, see formula lecture notes (third lecture by Julia Wiktor, on RDF formula)

counts = histogram[0]            #histogram[0] gives us an array with binnes values, notice histogram[1] gives the bin_edges, and you can also do hist[][] for indexing since its 2 dimensional.
counts = np.array(counts)     
counts = np.append(counts,1)  # i just use this to get the right dimensions when dividing in the rdf but we can probably solve it in a better way than this. 

dr = histogram[1][1]-histogram[1][0] #thickness of a shell, the indexing is used to extract bin_edges. Notice for example histogram[1][1] gives second bin_edge and [1][0] gives first bin_edge. 
print("thickness of a shell dr",dr)
r = histogram[1] - dr #this gives me bin_edge and subtracts the thickness from the shell so effectively we get the radius of each shell. imagine we have for example 4th shell. we subtract thicness to get radius of third shell and so on.


#r, dr and rho gives us all info we need to calculate volume of each shell. We now calculate the rdf by taking counts divided by the increasing shell volumes.

rdf = counts / ((4/3*np.pi*r**2*rho*dr)*len(atoms) ) #formula from third lecture by Julia Wiktor, on RDF)

integral = 0
for i, r_value in enumerate(r):
    if r_value <= 3.14:
        integral += rdf[i] * dr 
print("Integral value up to r = 3.14:", integral)

plt.plot(r,rdf)
plt.show()





#----------------------below is not important----------------------------#


#print("RDF",rdf)
#print("Volume",V)
#print("number of particles",N)
#print("Histogram[0] and also counts",counts)

#print("r is the radius of a shell: ",r)
#hist[0] gives an array with binned values
#hist[1] gives an array with bin_edges.
#with [0][1] we get binned values, numbers in first bin. 
#print(np.shape(r))
#print(np.shape(counts))
