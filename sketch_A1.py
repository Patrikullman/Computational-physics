import numpy as np

alpha_1 = 0.297104
alpha_2 = 1.236745
alpha_3 = 5.749982
alpha_4 = 38.216677

alpha_vec = [alpha_1, alpha_2, alpha_3, alpha_4]
alpha = np.zeros((4,1))
alpha = alpha_vec
Cp_vec = [1, 1, 1, 1]
Cp = np.zeros((4,1))
Cp = Cp_vec
norm = np.linalg.norm(Cp)
Cp = Cp / norm # Normalizing start guess


S = np.zeros((4,4))
T = np.zeros((4,4))
V = np.zeros((4,4))
h = np.zeros((4,4))
Q = np.zeros((4,4,4,4))
F = np.zeros((4,4))

for p in range(4):
    for q in range(4):
        S[p,q] = ( np.pi / (alpha[p] + alpha[q]) ) ** (3/2)
        #print(S[p,q])
        
for p in range(4):
    for q in range(4):
        T[p,q] = ( 3 * alpha[p]*alpha[q]* (np.pi**(3/2)) ) / ((alpha[p] + alpha[q])**(5/2))
        #print(T[p,q])

for p in range(4):
    for q in range(4):
        4*np.pi / (alpha[p] + alpha[q]) # Double check this integral!!!


h = T + V 
#print(h)

for p in range(4):
    for q in range(4):
        for r in range(4):
            for s in range(4):
                num = 2*(np.pi**(5/2)) # Numerator
                den = ((alpha[p] + alpha[q])*(alpha[r] + alpha[s])) * ((alpha[p] + alpha[q] + alpha[r] + alpha[s])**(1/2)) #Denominator
                Q[p,r,q,s] = num / den


for p in range(4):
    for q in range(4):
        term = np.zeros((4,4))
        for r in range(4):
            for s in range(4):
                term[p,q] += Q[p,q,r,s]*Cp[r]*Cp[s]
        F[p,q] = h[p,q] + term[p,q]

print(F)