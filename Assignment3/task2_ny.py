from ase.io import read
from ase import Atoms
from ase.visualize import view
from gpaw import GPAW
from ase.md.npt import NPT
from ase.io.trajectory import TrajectoryReader
import numpy as np
import matplotlib.pyplot as plt


traj = TrajectoryReader('NaCluster24.traj') 
all_distances = []
for i in range(13999):
    #print(traj[i][72]) # First index is timestep, second index is atom
    atoms_current_time_frame = traj[i]#[72]
    #print(atoms_current_time_frame[i])
    current_pos = atoms_current_time_frame.get_positions()
    Na_pos = current_pos[72]
    #print(Na_pos)

    O_pos_array = []
    for j in range(24):
        O_pos = current_pos[j]
        O_pos_array.append(O_pos)

    O_pos_array = np.array(O_pos_array)
    #print(O_pos_array.shape)
    L = 8.95
    #dist_array = []
    for j in range(24):
        diff_x = (Na_pos[0] - O_pos_array[j][0]) - L * np.round((Na_pos[0] - O_pos_array[j][0])/L) #calculates all distances, and also takes minimum image convention into consideration.
        diff_y = (Na_pos[1] - O_pos_array[j][1]) - L * np.round((Na_pos[1] - O_pos_array[j][1])/L)
        diff_z = (Na_pos[2] - O_pos_array[j][2]) - L * np.round((Na_pos[2] - O_pos_array[j][2])/L)
        dist = np.sqrt(diff_x**2 + diff_y**2 + diff_z**2)
        #dist_array.append(dist)
        all_distances.append(dist)
    #dist_array = np.array(dist_array)
all_distances = np.array(all_distances)
    
print(all_distances.shape)

bins = 1000
N = 24

counts, bin_edges = np.histogram(all_distances, bins)
bin_centers = (bin_edges[:-1] + bin_edges[1:])/2

rdf = counts/(4/3 * np.pi * (bin_centers**3))
plt.plot(bin_centers, rdf)
plt.show()


x = []
y = []
i = 0
while bin_centers[i] < 3.14:
    x.append(bin_centers[i])
    y.append(rdf[i])
    i += 1
    

def integrate(x, y):
    area = np.trapz(y=y, x=x)
    return area

area = integrate(x,y)
print(area)



# some comments :   we only need to calculcate RDF for the relevant species.  that means we only need RDF for sodium-oxygen. We dont need to consider hydrogen. 



