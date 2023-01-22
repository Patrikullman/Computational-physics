import numpy as np
from numpy.linalg import solve
import matplotlib.pyplot as plt

# Define the grid
N = 101
r_min = 0
r_max = 5
r = np.linspace(r_min, r_max, N)
dr = r[1] - r[0]
print(dr)


A = np.zeros((N,N))
b = np.zeros((N,1))


# Construct array A for finite difference approximation
for i in range(N):
    for j in range(N):
        if i == j:
            A[i,j] = -2
        elif (j == i-1):
            A[i,j] = 1
        elif (j == i+1):
            A[i,j] = 1
        else:
            A[i,j] = 0

A[0,1] = 0
A[N-1, N-2] = 0

# Define right hand side of equation
def ns(r):
    return (np.exp(-2*r) * (r**2)) / (np.pi)

def u(r):
    return np.sqrt(4*np.pi*ns(r))*r

rho = - (u(r)**2)/r * (dr**2)
b = rho

# Implement boundary condition
A[0,0] = 1
A[N-1,N-1] = 1
b[0] = 0
b[N-1] = 0

print(A)
print(b)

# Solve the system of equations
U = solve(A,b)
print(U)

plt.figure(1)
plt.plot(r, U)
plt.show()

# Plotting the Hartree potential to check the result of our program

# Define the Hartree potential
def VH(r):
    return 1/r - (1 + 1/r)*np.exp(-2*r)

V = VH(r)

# Plot the Hartree potential
plt.figure(2)
plt.plot(r, V)
plt.show()