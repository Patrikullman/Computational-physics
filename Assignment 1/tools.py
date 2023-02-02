import numpy as np
from scipy.linalg import eigh 
d = 4

def Q_fun(alpha):
    Q_array = np.zeros((d, d, d, d))
    for p in range(d):
        for r in range(d):
            for q in range(d):
                for s in range(d):
                    Q_array[p,r,q,s] = 2*np.pi**(5/2)/((alpha[p] +\
                    alpha[q])*(alpha[r] + alpha[s])*np.sqrt(alpha[p]+\
                    alpha[q]+alpha[r]+alpha[s]))
    return Q_array

def h_fun(alpha):
    h_array = np.zeros((d, d))
    for p in range(d):
        for q in range(d):
            h_array[p,q] = 3*alpha[p]*alpha[q]*np.pi**(3/2)/(alpha[p] + alpha[q])**(5/2)\
            - 4*np.pi/(alpha[p] + alpha[q]) 
    return h_array

def S_fun(alpha):
    S_array = np.zeros((d, d))
    for p in range(d):
        for q in range(d):
            S_array[p,q] = (np.pi/(alpha[p] + alpha[q]))**(3/2)
    return S_array

def F_fun(h, Q, C):
    F_array = np.zeros((d, d))
    for p in range(d):
        for q in range(d):
            term = 0
            for r in range(d):
                for s in range(d):
                    term += Q[p, r, q, s]*C[r]*C[s]
            F_array[p, q] = h[p, q] + term
    return F_array

def E_G_fun(C, h, Q):
    E_G_temp = 0
    for p in range(d):
        for q in range(d):
            E_G_temp += 2*C[p]*C[q]*h[p, q]
            for r in range(d):
                for s in range(d):
                    E_G_temp += Q[p,r,q,s]*C[p]*C[q]*C[r]*C[s]
    return E_G_temp

def normalize_C(C, S):
    norm_term = 0 
    for p in range(d):
        for q in range(d):
            norm_term += C[p]*S[p, q]*C[q]
    return C/np.sqrt(norm_term)

def psi_H_1s(r):
    return (1/np.sqrt(np.pi))*np.exp(-r)

def V_H(r):
    return 1/r - (1 + 1/r)*np.exp(-2*r)


def Kohn_Sham(V_tot, a, b, d, Z):
    h = (b-a)/d
    r = np.zeros((d-2))
    for i in range(0, d-2):
        r[i] = a + (i + 1)*h
        
    # Define the matrix A
    A = np.zeros((d-2, d-2))
    for i in range(d-2):
        for j in range(d-2):
            if(i == j):
                A[i, j] = 1/(h**2) - Z/(r[i]) + V_tot[i]
            elif(i == j - 1 or i == j + 1):
                A[i, j] = -1/(2*h**2)
            else:
                A[i, j] = 0

    # Solve the eigenvalue problem
    eig_vals, eig_vecs = eigh(A)

    # eigenvalues and eigenvectors
    E = eig_vals[0]
    E_vec = eig_vecs[:, 0]

    # Wavefunction
    phi = E_vec/(r*np.sqrt(4*np.pi))
    # Normalize wavefunction
    term =np.trapz(phi**2*4*np.pi*r**2,x = r, dx = h)
    phi = phi/np.sqrt(term)
    
    return E, phi

# 
def get_V_sH(phi_guess, a, b, d):
    h = (b-a)/d
    r = np.zeros((d-2))
    for i in range(0, d-2):
        r[i] = a + (i + 1)*h

    u = np.sqrt(4*np.pi)*r*phi_guess

    # Ax=B
    A = np.zeros((d-2, d-2))
    B = np.zeros((d-2))
    for i in range(d-2):
        B[i] = h**2*u[i]*u[i]/r[i]
        for j in range(d-2):
            if(i == j):
                A[i, j] = 2
            elif(i == j - 1 or i == j + 1):
                A[i, j] = -1
            else:
                A[i, j] = 0

    # Solve the matrix equation
    U = np.linalg.solve(A, B)

    V_sH = U/(r) + 1/b
    return V_sH

# Calculate V_xc from eq. (22-27)
def get_V_xc(phi):
    A = 0.0311; B = -0.048; C = 0.002; D = -0.0116
    gamma = -0.1423; beta_1 = 1.0529; beta_2 = 0.3334

    n = 2*np.abs(phi)**2
    r_s = (3/(4*np.pi*n))**(1/3)
    eps_x = -(3/4*(3*n/np.pi)**(1/3))
    eps_x_der = -3/(4*np.pi)*(3*n/np.pi)**(-2/3)
    eps_c = np.zeros((len(n)))
    eps_c_der = np.zeros((len(n)))
    for i in range(len(n)):
        if(r_s[i] < 1):
            eps_c[i] = A*np.log(r_s[i]) + B + \
                C*r_s[i]*np.log(r_s[i]) + D*r_s[i]
            eps_c_der[i] = -(4**2*np.pi**2/(3**5))**(1/3)*r_s[i]*(A/r_s[i]-\
                 C*(np.log(r_s[i]) + 1) + D)
        else:
            eps_c[i] = gamma / (1 + beta_1*np.sqrt(r_s[i]) + beta_2*r_s[i])
            eps_c_der[i] = (4**2*np.pi**2/(3**5))**(1/3)*r_s[i]*(gamma*((beta_2+beta_1/\
                (2*r_s[i]**(1/2))))/((1+beta_1*r_s[i]**(1/2)+beta_2*r_s[i])**2))

    eps_xc = eps_x + eps_c
    eps_xc_der = eps_x_der + eps_c_der

    V_xc = eps_xc + n*eps_xc_der
    return V_xc, eps_xc