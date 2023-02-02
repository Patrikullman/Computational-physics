import numpy as np
import matplotlib.pyplot as plt
from tools import*

# Initialize the r-vector constrained with a constant point density
density = 40; a = 0; b = 12
d = (b-a)*density
h = 1/density
Z = 2
r = np.zeros((d-2))
for i in range(0, d-2):
    r[i] = a + (i + 1)*h

E_0 = 0; E_old =0; E_del = 1; eps = 0
phi_guess = np.ones((d-2))
while(E_del > 1e-5):
    # Get V_sH by solving EVP
    V_sH = get_V_sH(phi_guess, a, b, d)
    V_H = 2*V_sH

    # Get V_xc with only exchange term 
    n = 2*np.abs(phi_guess)**2
    V_x = -(3/4*(3*n/np.pi)**(1/3)) - 3*n/(4*np.pi)*(3*n/np.pi)**(-2/3)
    V_tot = V_H + V_x

    # Solve Kohn-Sham EVP
    eps, phi_guess = Kohn_Sham(V_tot, a, b, d, Z)

    # Calculate ground state energy
    u = phi_guess*r*np.sqrt(4*np.pi)
    term = np.trapz(u**2*((1/2)*V_H + V_x + (3/4*(3*n/np.pi)**(1/3))), x = r, dx = h)
    E_0 = 2*eps - 2*term

    E_del = np.abs(E_old - E_0)
    E_old = E_0
print("eps: " + str(eps))
print("E_0: " + str(E_0))

term =np.trapz(phi_guess**2*4*np.pi*r**2,x = r, dx = h)
phi_guess = phi_guess/np.sqrt(term)

fig, ax = plt.subplots()

phi_helium = np.genfromtxt('helium.csv', delimiter=',', skip_header=0)

ax.set_ylabel(r'Wavefunction $\psi(r)$', fontsize=15)
ax.set_xlabel(r'Radius $r$ (a.u.)', fontsize=15)
ax.grid()
ax.plot(r, phi_guess, 'b', linewidth = 4, label = r"$\psi(r)$ using DFT with $V_x$")
ax.plot(r, phi_helium, '--r', linewidth = 4, label = r"$\psi(r)$ from Task 1")
ax.tick_params(axis = 'both', which = 'major', labelsize = 12)
fig.tight_layout()
ax.legend(prop={'size': 17})

fig.savefig("a1-5.pdf")