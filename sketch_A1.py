import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt

# Set the four alpha values
alpha_1 = 0.297104
alpha_2 = 1.236745
alpha_3 = 5.749982
alpha_4 = 38.216677

# Make alpha vector
alpha_vec = [alpha_1, alpha_2, alpha_3, alpha_4]
alpha = np.zeros((4,1))
alpha = alpha_vec

# Make Cp vector with start guesses
Cp_vec = [1, 1, 1, 1]
Cp = np.zeros((4,1))
Cp = Cp_vec

# Define some matrices filled with zeros that will be filled with values later
S = np.zeros((4,4))
T = np.zeros((4,4))
V = np.zeros((4,4))
h = np.zeros((4,4))
Q = np.zeros((4,4,4,4))
F = np.zeros((4,4))

# Insert values in S matrix
for p in range(4):
    for q in range(4):
        S[p,q] = ( np.pi / (alpha[p] + alpha[q]) ) ** (3/2)
        
# S matrix defined, now Cp can be normalized
# Normalize Cp accordin to the condition sum_ij Cp[p]*S[p,q]*Cp[q] = 1
Cp = Cp / np.sqrt(np.dot(Cp, np.dot(S,Cp)))

# Insert values in the T matrix (holding kinetic energy values)
for p in range(4):
    for q in range(4):
        T[p,q] = ( 3 * alpha[p]*alpha[q]* (np.pi**(3/2)) ) / ((alpha[p] + alpha[q])**(5/2))

# Insert values in the V matrix (holding potential energy values)
for p in range(4):
    for q in range(4):
        V[p,q] = -4*np.pi / (alpha[p] + alpha[q]) # Double check this integral!!!

# Sum of T and V gives the total Hamiltonian h
h = T + V 

# Insert values in the 4D Q matrix
for p in range(4):
    for q in range(4):
        for r in range(4):
            for s in range(4):
                num = 2*(np.pi**(5/2)) # Numerator
                den = ((alpha[p] + alpha[q])*(alpha[r] + alpha[s])) * ((alpha[p] + alpha[q] + alpha[r] + alpha[s])**(1/2)) #Denominator
                Q[p,r,q,s] = num / den

# Insert values in the F (Fock) matrix
for p in range(4):
    for q in range(4):
        term = np.zeros((4,4))
        for r in range(4):
            for s in range(4):
                term[p,q] += Q[p,q,r,s]*Cp[r]*Cp[s]
        F[p,q] = h[p,q] + term[p,q]

i = 0
energy = []
ite = []
while i < 100:
    eigenvalues, eigenvectors = eigh(F, S)
    #print("Eigenvalues:", eigenvalues)
    #print("Eigenvectors:", eigenvectors)

    min_eig_index = np.argmin(eigenvalues) #finding index of lowest eigenvalue because we need to extract corresponding eigenvector.
    #print(min_eig_index)
    
    
    min_eig_vec = eigenvectors[min_eig_index,:] #this eigenvector corresponds to lowest eigenvalue, that is the ground state and our "new Cp" that is used to update F matrix.
    #print(min_eig_vec)
    Cp = min_eig_vec 
    Cp = Cp / np.sqrt(np.dot(Cp, np.dot(S,Cp)))

    sum1 = 0
    sum2 = 0

    for p in range(4):
        for q in range(4):
            sum1 += Cp[p]*Cp[q]*h[p,q]

    for p in range(4):
        for q in range(4):
            for r in range(4):
                for s in range(4):
                    sum2 += Q[p,q,r,s]*Cp[p]*Cp[q]*Cp[r]*Cp[s]
    
    # Calculating ground state energy
    Eg = 2*sum1 + sum2
    print(Eg)
    energy.append(Eg)
    ite.append(i)

    # Update the F matrix with the new normalized Cp valuess
    for p in range(4):
        for q in range(4):
            term = np.zeros((4,4))
            for r in range(4):
                for s in range(4):
                    term[p,q] += Q[p,q,r,s]*Cp[r]*Cp[s]
            F[p,q] = h[p,q] + term[p,q]

    i += 1

plt.plot(ite, energy)
plt.show()

