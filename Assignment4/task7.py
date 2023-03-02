from ase.build import molecule
import numpy as np
from ase.thermochemistry import IdealGasThermo
import matplotlib.pyplot as plt




kB = 8.617333E-5 # Boltzmann (eV/K)
nu = 1E+12 # (s^-1)
T = np.linspace(100, 2000, 101)

######### This is for the Au surface ###########

#E_gas_CO = -14.210 # (eV)
E_adsorbate_CO = -0.137 # (eV)
E_adsorbate_O = -0.057 # (eV)
E_activation = 0.278 # (eV)

S_gas_O = np.zeros((len(T))) # Array holding entropy for different temperature
S_gas_CO = np.zeros((len(T)))

CO = molecule('CO', vacuum = 12)
O2 = molecule('O2', vacuum = 12)

thermoCO = IdealGasThermo(vib_energies = [0.22],         
                        potentialenergy = 0,
                        geometry = 'linear', atoms = CO,
                        symmetrynumber = 1, spin = 0)  #O2 :: spin= 1 , symmetry  = 2  #CO :: spin = 0 , symmetry =1 

thermoO2 = IdealGasThermo(vib_energies = [0.204],         
                        potentialenergy = 0,
                        geometry= 'linear', atoms = O2,
                        symmetrynumber = 2, spin = 1) 

for i in range(len(T)):    #we need entropies for each temperature step 
   S_gas_CO[i] = thermoCO.get_entropy(temperature = T[i], pressure = 101325, verbose = False)                #should we calculate the DFT and entropy for every temperature step ? 
   S_gas_O[i] = thermoO2.get_entropy(temperature = T[i], pressure = 101325, verbose = False)


Keq_O = np.zeros((len(T)))
Keq_CO = np.zeros((len(T)))


for i in range(len(T)):
    Keq_O[i] = np.exp(-S_gas_O[i] / kB) * np.exp((E_adsorbate_O) / (kB * T[i]))
    Keq_CO[i] += np.exp(-S_gas_CO[i] / kB) * np.exp((E_adsorbate_CO) / (kB * T[i]))
    
print(Keq_O)
print(Keq_CO)

theta_CO = np.zeros((len(T)))
theta_O = np.zeros((len(T)))
theta_CO2 = np.zeros((len(T)))

for i in range(len(T)):
    #theta_CO[i] = (Keq_CO[i]*101325) / (1 + Keq_CO[i]*101325)
    #theta_O[i] = (Keq_O[i]*101325) / (1 + Keq_O[i]*101325)

    theta_CO[i] = (Keq_CO[i]*101325) / (1 + Keq_CO[i]*101325 + Keq_O[i]*101325)
    theta_O[i] = (Keq_O[i]*101325) / (1 + Keq_CO[i]*101325 + Keq_O[i]*101325)

    #theta_CO2[i] = theta_O[i] * theta_CO[i]**2 # not sure if the ^2 is on the right molecule




r = np.zeros((len(T)))

for i in range(len(T)):
    r[i] = theta_CO2[i]*nu*np.exp(-E_activation / (kB*T[i]))

plt.figure(1)
plt.plot(T, theta_CO, 'b')
plt.title('theta CO on Au')

plt.figure(2)
plt.plot(T, theta_O, 'b')
plt.title('theta O on Au')



######### This is for the Pt surface ###########

E_adsorbate_CO = -1.659 # (eV)
E_adsorbate_O = -1.171 # (eV)
E_activation = 1.069 # (eV)

S_gas_O = np.zeros((len(T))) # Array holding entropy for different temperature
S_gas_CO = np.zeros((len(T)))

CO = molecule('CO', vacuum = 12)
O2 = molecule('O2', vacuum = 12)

thermoCO = IdealGasThermo(vib_energies = [0.22],         
                        potentialenergy = 0,
                        geometry = 'linear', atoms = CO,
                        symmetrynumber = 1, spin = 0)  #O2 :: spin= 1 , symmetry  = 2  #CO :: spin = 0 , symmetry =1 

thermoO2 = IdealGasThermo(vib_energies = [0.204],         
                        potentialenergy = 0,
                        geometry= 'linear', atoms = O2,
                        symmetrynumber = 2, spin = 1) 

for i in range(len(T)):    #we need entropies for each temperature step 
   S_gas_CO[i] = thermoCO.get_entropy(temperature = T[i], pressure = 101325, verbose = False)                #should we calculate the DFT and entropy for every temperature step ? 
   S_gas_O[i] = thermoO2.get_entropy(temperature = T[i], pressure = 101325, verbose = False)



print(S_gas_O)
print(S_gas_CO)

Keq_O = np.zeros((len(T)))
Keq_CO = np.zeros((len(T)))


for i in range(len(T)):
    Keq_O[i] = np.exp(-S_gas_O[i] / kB) * np.exp((E_adsorbate_O) / (kB * T[i]))
    Keq_CO[i] += np.exp(-S_gas_CO[i] / kB) * np.exp((E_adsorbate_CO) / (kB * T[i]))
    
print(Keq_O)
print(Keq_CO)

theta_CO = np.zeros((len(T)))
theta_O = np.zeros((len(T)))
theta_CO2 = np.zeros((len(T)))

for i in range(len(T)):
    #theta_CO[i] = (Keq_CO[i]*101325) / (1 + Keq_CO[i]*101325)
    #theta_O[i] = (Keq_O[i]*101325) / (1 + Keq_O[i]*101325)
    #theta_CO2[i] = theta_O[i]**2 * theta_CO[i] # not sure if the ^2 is on the right molecule
    theta_CO[i] = (Keq_CO[i]*101325) / (1 + Keq_CO[i]*101325 + Keq_O[i]*101325)
    theta_O[i] = (Keq_O[i]*101325) / (1 + Keq_CO[i]*101325 + Keq_O[i]*101325)

plt.figure(3)
plt.plot(T, theta_CO, 'r')
plt.title('theta CO on Pt')

plt.figure(4)
plt.plot(T, theta_O, 'r')
plt.title('theta O on Pt')


######### This is for the Rh surface ###########

E_adsorbate_CO = -1.75 # (eV)
E_adsorbate_O = -1.735 # (eV)
E_activation = 1.266 # (eV)


S_gas_O = np.zeros((len(T))) # Array holding entropy for different temperature
S_gas_CO = np.zeros((len(T)))

CO = molecule('CO', vacuum = 12)
O2 = molecule('O2', vacuum = 12)

thermoCO = IdealGasThermo(vib_energies = [0.22],         
                        potentialenergy = 0,
                        geometry = 'linear', atoms = CO,
                        symmetrynumber = 1, spin = 0)  #O2 :: spin= 1 , symmetry  = 2  #CO :: spin = 0 , symmetry =1 

thermoO2 = IdealGasThermo(vib_energies = [0.204],         
                        potentialenergy = 0,
                        geometry= 'linear', atoms = O2,
                        symmetrynumber = 2, spin = 1) 

for i in range(len(T)):    #we need entropies for each temperature step 
   S_gas_CO[i] = thermoCO.get_entropy(temperature = T[i], pressure = 101325, verbose = False)                #should we calculate the DFT and entropy for every temperature step ? 
   S_gas_O[i] = thermoO2.get_entropy(temperature = T[i], pressure = 101325, verbose = False)



print(S_gas_O)
print(S_gas_CO)

Keq_O = np.zeros((len(T)))
Keq_CO = np.zeros((len(T)))


for i in range(len(T)):
    Keq_O[i] = np.exp(-S_gas_O[i] / kB) * np.exp((E_adsorbate_O) / (kB * T[i]))
    Keq_CO[i] += np.exp(-S_gas_CO[i] / kB) * np.exp((E_adsorbate_CO) / (kB * T[i]))
    
print(Keq_O)
print(Keq_CO)

theta_CO = np.zeros((len(T)))
theta_O = np.zeros((len(T)))
theta_CO2 = np.zeros((len(T)))

for i in range(len(T)):
    #theta_CO[i] = (Keq_CO[i]*101325) / (1 + Keq_CO[i]*101325)
    #theta_O[i] = (Keq_O[i]*101325) / (1 + Keq_O[i]*101325)
    #theta_CO2[i] = theta_O[i]**2 * theta_CO[i] # not sure if the ^2 is on the right molecule
    theta_CO[i] = (Keq_CO[i]*101325) / (1 + Keq_CO[i]*101325 + Keq_O[i]*101325)
    theta_O[i] = (Keq_O[i]*101325) / (1 + Keq_CO[i]*101325 + Keq_O[i]*101325)



plt.figure(5)
plt.plot(T, theta_CO, 'lime')
plt.title('theta CO on Rh')

plt.figure(6)
plt.plot(T, theta_O, 'lime')
plt.title('theta O on Rh')

plt.show()


