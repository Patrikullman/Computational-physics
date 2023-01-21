import numpy as np
from scipy.linalg import eig
import numpy as np
alpha_1 = 0.297104 #defining alpha values
alpha_2 = 1.236745
alpha_3 = 5.749982
alpha_4 = 38.216677

alpha_vec = [alpha_1, alpha_2, alpha_3, alpha_4] #putting alpha values in an array
alpha = np.zeros((4, 1))
alpha = alpha_vec
Cp_vec = [1, 1, 1, 1]
Cp = np.zeros((4, 1))
Cp = Cp_vec
norm = np.linalg.norm(Cp)
Cp = Cp / norm  # Normalizing start guess

S = np.zeros((4, 4))
T = np.zeros((4, 4))
V = np.zeros((4, 4))
h = np.zeros((4, 4))
Q = np.zeros((4, 4, 4, 4))
F = np.zeros((4, 4))

for p in range(4):
    for q in range(4):
        S[p, q] = (np.pi / (alpha[p] + alpha[q])) ** (3 / 2)  #constructing S matrix
        # print(S[p,q])

for p in range(4):
    for q in range(4):
        T[p, q] = (3 * alpha[p] * alpha[q] * (np.pi ** (3 / 2))) / ((alpha[p] + alpha[q]) ** (5 / 2))
        # print(T[p,q])

for p in range(4):
    for q in range(4):
       V = 4 * np.pi / (alpha[p] + alpha[q])  # Double check this integral!!!

h = T + V
# print(h)

for p in range(4):
    for q in range(4):
        for r in range(4):
            for s in range(4):
                num = 2 * (np.pi ** (5 / 2))  # Numerator
                den = ((alpha[p] + alpha[q]) * (alpha[r] + alpha[s])) * (
                            (alpha[p] + alpha[q] + alpha[r] + alpha[s]) ** (1 / 2))  # Denominator
                Q[p, r, q, s] = num / den

for p in range(4):
    for q in range(4):
        term = np.zeros((4, 4))
        for r in range(4):
            for s in range(4):
                term[p, q] += Q[p, q, r, s] * Cp[r] * Cp[s]
        F[p, q] = h[p, q] + term[p, q]



#we know that Cp will change, it will converge to some values? somehow C_r and C_s are updated and put in the new Cp vector iteratively.
#vi behöver en while loop som jämför E värden och när skilladen är 10^-5 så hoppar den ut och printar ut E.
#vi behöver en Cp_new och en Cp_old eller liknande för Cp kommer ju att uppdateras iterativt.


#take the minimum eigenvalue that we get from np.linalg.eig(F,S), this is the eigenvalue corresponding to the C eigenvector we want. This new
#C vector goes into the F matrix, and a new eigenvalue and eigenvectors are computed. This continues until we reach a eigenvalue close to the energy that we are looking for

#It is the F matrix that changes every time, because it contains the values for Cr and Cs and those change after every iteration.

import matplotlib.pyplot as plt

#-------------------------------------------------------------------------
i = 1
while i < 15: #E_new-E_old == 10**(-5):
    for p in range(4):
        for q in range(4):
            term = np.zeros((4, 4))
            for r in range(4):
                for s in range(4):
                    term[p, q] += Q[p, q, r, s] * Cp[r] * Cp[s]
            F[p, q] = h[p, q] + term[p, q]
    #print(F)

    eigenvalues, eigenvectors = eig(F, S) #solves generalized eigenvalue problem of FC = E'SC, we obtain four eigenvalues and four eigenvectors.
    #print("Eigenvalues:", eigenvalues)
    #print("Eigenvectors:", eigenvectors)

    min_eig_index = np.argmin(eigenvalues) #finding index of lowest eigenvalue because we need to extract corresponding eigenvector.
    #print(min_eig_index)
    min_eig_vec = eigenvectors[min_eig_index,:] #this eigenvector corresponds to lowest eigenvalue, that is the ground state and our "new Cp" that is used to update F matrix.
    Cp = min_eig_vec #setting our Cp guess to our new solution.
    #print("this is Cp",Cp)
    Cp = Cp / np.linalg.norm(Cp) #we normalize our new Cp after each iteration.
    print("this is Cp normalized:",Cp)

    i += 1 #we loop up until C is stable, this happens after 10 iterations... Then we can use our C to calculate the actual ground state energy i think.
    #and the wavefunction is then given becuase we calculated our weights C so just plug them into equation 4.12 (page 49).








