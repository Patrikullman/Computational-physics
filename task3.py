import numpy as np
from scipy.linalg import eigh

# Define the grid
N = 101
r_min = 0
r_max = 5
r = np.linspace(r_min, r_max, N)
dr = r[1] - r[0]

# Create a Hamiltonian matrix
H = np.zeros((N,N)) 

# Define the Hartree potential
def V_H(r):
    return 1/r # Hartree potential for Hydrogen

# Define the exchange potential
def V_x(r):
    return 0 # Let this be zero at the moment

# Define the correlation potentihal
def V_c(r):
    return 0 # Let this be zero at the moment


# Fill the Hamiltonian matrix based on the finite differences
for i in range(1,N-1):
    H[i, i-1] = -1 / (2 * dr**2)
    H[i, i] = 1 / dr**2 + V_H(r[i]) + V_x(r[i]) + V_c(r[i])
    H[i, i+1] = -1 / (2 * dr**2)

# Implement the boundary conditions
H[0,0] = 1
H[-1,-1] = 1

#print(H)

# Solve the eigenvalue problem
E, C = eigh(H)
print(E)










