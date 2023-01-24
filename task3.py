import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt

# Define the grid
N = 1000
r_min = 0
r_max = 10
r = np.zeros((N))
dr = (r_max-r_min)/N
for i in range(N):
    r[i] = r_min + i*dr
#r = np.linspace(r_min, r_max, N)
#dr = r[1] - r[0]



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
    H[i, i] = 1 / dr**2 -1/r[i]#- 2/r[i] + V_H(r[i]) + V_x(r[i]) + V_c(r[i])
    H[i, i+1] = -1 / (2 * dr**2)

# Implement the boundary conditions
H[0,0] = 0
H[-1,-1] = 0
H[1,0] = 0


#print(H)

#print(H)

# Solve the eigenvalue problem
E, C = eigh(H)
print(E[0]*27.21132)
print(E[0])

# Solving for the wavefunction phi
u = C[:,0]
#print(u)
phi = u / (np.sqrt(4*np.pi)*r)
phi = np.delete(phi, 0)
r = np.delete(r, 0)
#print(phi)


integral = np.trapz(4*np.pi*r**2*phi**2, r)
#print(integral)
#summa = np.sum(phi)
phi = phi / np.sqrt(integral)
#print(sum(phi))


# The normalized hydrogen ground state wavefunction
def psi(r):
    return 2*np.exp(-r)/(2*np.sqrt(np.pi))    

psi_hydrogen = psi(r)
plt.plot(r, phi, label = 'Solution of Kohn-Sham')
plt.plot(r, psi_hydrogen, 'r--', label = 'Wavefunction of hydrogen ground state')
plt.xlabel('Radius (atomic units)', fontsize = 12)
plt.ylabel('$\phi$', fontsize = 12)
plt.legend(loc = 1, fontsize = 12)
plt.show()










