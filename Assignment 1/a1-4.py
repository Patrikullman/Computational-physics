import numpy as np
from scipy.linalg import eigh 
import matplotlib.pyplot as plt
from tools import*

d = 1000; a = 0; b = 10
h = (b-a)/d
Z = 2
r = np.zeros((d-2))
for i in range(0, d-2):
    r[i] = a + (i + 1)*h

print(r)

E_0 = 0; E_old =0; E_del = 1
phi_guess = np.ones((d-2))
while(E_del > 1e-5):

    V_sH = get_V_sH(phi_guess, a, b, d)

    eps, phi_guess = Kohn_Sham(V_sH, a, b, d, Z)

    term = np.trapz((1/2)*V_sH*phi_guess**2*r**2*4*np.pi, x = r, dx = h)

    E_0 = 2*eps - 2*term

    E_del = np.abs(E_old - E_0)
    E_old = E_0
    print(E_0)
print("-----")

'''
density = 40
a = 0
b = 5
d = (b-a)*density
h = 1/density
Z = 2
r = np.zeros((d-2))
for i in range(0, d-2):
    r[i] = a + (i + 1)*h

phi_guess = np.ones((d-2))

n_b = 40
b_array = np.zeros((n_b))
E_0_array = np.zeros((n_b))

for i in range(n_b):
    b = 1 + i*0.2
    d = int((b-a)*density)
    h = 1/density
    r = np.zeros((d-2))
    for j in range(0, d-2):
        r[j] = a + (j + 1)*h
    phi_guess = np.ones((d-2))

    E_0 = 0
    E_old =0
    E_del = 1
    while(E_del > 1e-5):
        b_array[i] = b

        V_sH = get_V_sH(phi_guess, a, b, d)

        vec = Kohn_Sham(V_sH, a, b, d, Z)

        eps = vec[0]
        phi_guess = vec[1]

        term = np.trapz((1/2)*V_sH*phi_guess**2*r**2*4*np.pi, x = r, dx = h)

        E_0 = 2*eps - 2*term

        E_del = np.abs(E_old - E_0)
        E_old = E_0
    E_0_array[i] = E_0
    print("b: " + str(b) + ": " + str(E_0))

fig, ax = plt.subplots()

ax.set_ylabel(r'Ground state energy (a.u.)', fontsize=18)
ax.set_xlabel(r'$r_{max}$ (a.u.)', fontsize=18)
ax.grid()
ax.plot(b_array, E_0_array, 'b', linewidth = 3, label = r"Ground state energy $E_0$")
ax.tick_params(axis = 'both', which = 'major', labelsize = 12)
fig.tight_layout()
ax.legend(prop={'size': 17})

fig.savefig("a1-4_b.pdf")


density = 20
a = 0
b = 15
d = (b-a)*density
h = 1/density
Z = 2
r = np.zeros((d-2))
for i in range(0, d-2):
    r[i] = a + (i + 1)*h

phi_guess = np.ones((d-2))

n_d = 20
density_array = np.zeros((n_d))
E_0_array = np.zeros((n_d))

for i in range(n_d):
    density = 2 + i*2
    d = (b-a)*density
    h = 1/density
    r = np.zeros((d-2))
    for j in range(0, d-2):
        r[j] = a + (j + 1)*h
    phi_guess = np.ones((d-2))

    E_0 = 0
    E_old =0
    E_del = 1
    while(E_del > 1e-5):
        density_array[i] = density

        V_sH = get_V_sH(phi_guess, a, b, d)

        vec = Kohn_Sham(V_sH, a, b, d, Z)

        eps = vec[0]
        phi_guess = vec[1]

        term = np.trapz((1/2)*V_sH*phi_guess**2*r**2*4*np.pi, x = r, dx = h)

        E_0 = 2*eps - 2*term

        E_del = np.abs(E_old - E_0)
        E_old = E_0
    E_0_array[i] = E_0
    print("density: " + str(density) + ": " + str(E_0))

fig, ax = plt.subplots()

ax.set_ylabel(r'Ground state energy (a.u.)', fontsize=18)
ax.set_xlabel(r'Grid density $d$', fontsize=18)
ax.grid()
ax.plot(density_array, E_0_array, 'b', linewidth = 3, label = r"Ground state energy $E_0$")
ax.tick_params(axis = 'both', which = 'major', labelsize = 12)
fig.tight_layout()
ax.legend(prop={'size': 17})

fig.savefig("a1-4_d.pdf")

fig, ax = plt.subplots()

phi_helium = np.genfromtxt('helium.csv', delimiter=',', skip_header=0)

ax.set_ylabel(r'Wavefunction $\psi(r)$', fontsize=15)
ax.set_xlabel(r'Radius $r$ (a.u.)', fontsize=15)
ax.grid()
ax.plot(r, phi_guess, 'b', linewidth = 4, label = r"$\psi(r)$ using DFT with $V_{xc} = 0$")
ax.plot(r, phi_helium, '--r', linewidth = 4, label = r"$\psi(r)$ from Task 1")
ax.tick_params(axis = 'both', which = 'major', labelsize = 12)
fig.tight_layout()
ax.legend(prop={'size': 17})

fig.savefig("a1-4.pdf")
'''