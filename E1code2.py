#!/usr/bin/env python
###############################################################################
# E1code2
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

# skip_header skips the first
# row in data.csv
array = np.genfromtxt('powerspectrum.csv', delimiter=',', skip_header=1)

fig, ax = plt.subplots()
ax.plot(array[:, 1], array[:, 0])   #byt ut 0 till 1 och 1 till 0, och byt till powerspectrum.csv för att plotta powerspectrum

ax.set_xlabel('time (arb.unit)')
ax.set_ylabel('signal (arb.unit)')
ax.grid()

fig.savefig('signal.pdf')
plt.show()
