import numpy as np
from numpy.linalg import norm
from numpy.linalg import solve
from scipy.linalg import eigh
import matplotlib.pyplot as plt
from tools import *

# First we define some functions to get V_sH, the wave function phi and corresponding energy

# Define function to get V_sH (basically what task 2 does)
def get_V_sH(r, phi):
    N = len(phi)
    A = np.zeros((N,N))
    b = np.zeros((N,1))


    # Construct array A for finite difference approximation
    for i in range(N):
        for j in range(N):
            if (i == j):
                A[i,j] = -2
            elif (j == i-1):
                A[i,j] = 1
            elif (j == i+1):
                A[i,j] = 1
            else:
                A[i,j] = 0

    A[0,0] = 1
    A[0,1] = 0
    A[N-1, N-2] = 0
    A[N-1,N-1] = 1

    # Define right hand side of equation
    u = np.zeros((N))
    rho = np.zeros((N))
    for i in range(N):
        u[i] = np.sqrt(4*np.pi) * r[i] * phi[i]

    for i in range(N):
        rho[i] = - (u[i]**2)/r[i] * (dr**2)
    b = rho

    # Implement boundary condition
    A[0,0] = 1
    A[N-1,N-1] = 1
    b[0] = 0
    b[N-1] = 0

    U = solve(A,b)
    #print(b.shape)
    #r = np.delete(r, 0)
    U = U + r/r_max
    V_sH = U / r

    return V_sH


# Define function to get the wave function phi and the corresponding energy
def get_phi(r, V_H, V_x, V_c):
    N = len(V_H)

    # Create a Hamiltonian matrix
    H = np.zeros((N,N)) 


    # Fill the Hamiltonian matrix based on the finite differences
    # The BC is applied by letting the first and last row be filled with zeros
    for i in range(1,N-1):
        H[i, i-1] = -1 / (2 * dr**2)
        H[i, i] = 1 / dr**2 - 2/r[i] + V_H[i] + V_x[i] + V_c[i]
        H[i, i+1] = -1 / (2 * dr**2)
    H[1,0] = 0

    # Solve the eigenvalue problem
    E, C = eigh(H)
    #print(E[0]) # ground state energy in atomic units
    #print(E[0]*27.21132) # ground state energy in eV

    # Solving for the wavefunction phi
    u = C[:,0]
    phi = u / (np.sqrt(4*np.pi)*r)

    # Remove first element in phi and r corresponding to the division by zero term
    #phi = np.delete(phi, 0)
    #r = np.delete(r, 0)


    # Compute integral of phi**2 in spherical coordinates with the trapezoidal rule
    integral = np.trapz(4*np.pi*r**2*phi**2, r)

    # Normalize phi
    phi = phi / np.sqrt(integral)

    return E[0], phi


# Below is the program that solves the problem and calls the function

# Define grid
N = 1000
r_min = 0
r_max = 10
r = np.linspace(r_min, r_max, N)
dr = r[1] - r[0]

# Remove the r = 0 element to avoid division by zero in the functions
r = np.delete(r, 0)

# Start guess for wave function and normalize it
phi_start = np.ones((N-1))
phi_start = phi_start / norm(phi_start)
phi = phi_start
#print(phi.shape)

n = 3 / (4*np.pi*r**3)
epsilon_x = (-3 / 4) * ((3*n/np.pi)**(1/3))
V_x = epsilon_x + n*np.gradient(epsilon_x)




E_diff = 1
E_old = 1
while np.abs(E_diff) > 10**(-5):
    V_H_vector = 2*get_V_sH(r, phi)
    E, phi = get_phi(r, V_H_vector, V_x, np.zeros((N-1)))
    E_diff = (E - E_old)*27.21132 # Energy difference between new and old energy in eV
    E_old = E
    print(E_diff)
print(E)

plt.plot(r, phi)
plt.title('Ground state wave function for the He-atom', fontsize = 14)
plt.xlabel('Radius (atomic units)', fontsize = 12)
plt.ylabel(r'Wave function $\varphi (r)$', fontsize = 12)
plt.show()







