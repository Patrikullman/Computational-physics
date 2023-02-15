import numpy as np
import matplotlib.pyplot as plt
from ase.io import read

# Define number of iterations and bins
start_index = 10000
n_iterations = 2000
n_bins = 500

# Our histogram 
rdf = np.zeros((n_bins))
atoms = read ('cluster24.traj', index = '10000:')

for k in range (0, 24):
    for a in atoms :
        # Get distances for all the oxygen atoms
        get_dist_array = a.get_distances (k , indices = range (0 ,24) , mic = True )
        # Calculate histogram for that
        hist_temp = np.histogram (get_dist_array, n_bins, range = (1, 6))[0]
        # Add to our histogram
        rdf += hist_temp

# Normalise the histogram .
r = np.linspace (1, 6, n_bins)
dr = 6/n_bins
V = 8.956369831609006**3
N = 24
rho = N/V

for i in range(n_bins):
    rdf[i] = rdf[i]/(4*np.pi*r[i]**2*dr*rho*N)


# Divide with number of iterations, and save to csv file .
rdf = rdf/n_iterations


integral_upper_limit = 3.13
integral = 0
for i, r_value in enumerate(r):
    if r_value <= integral_upper_limit:
        integral += rdf[i] * r[i]**2 * dr
cord_num = 4*np.pi*rho*integral 
integral_up_lim_str = str(integral_upper_limit)
print('Integral value up to r = ' + integral_up_lim_str + ': ', cord_num)

plt.plot (r, rdf)
plt.fill_between(r, 0, rdf, where = (r <= 3.14),alpha = 0.5, color = 'grey', label = 'First solvation shell')
plt.xlabel('Radius (r)', fontsize = 12)
plt.ylabel('Counts ', fontsize = 12)
plt.title('cluster24.traj', fontsize = 14)
plt.legend(loc = 'upper right', fontsize = 12)
plt.show()


