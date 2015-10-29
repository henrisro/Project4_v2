# Makefile program for running Project4_d.cpp

import os

os.system('g++ -c Project4_d.cpp lib.cpp')
os.system('g++ -o Project4_d.out Project4_d.o lib.o -O2')
# Output file should be provided here as well:
os.system('./Project4_d.out Ising2D_d.txt')

# Read data from 'Ising2D_d.txt' and plot energy histogram,
# i.e. energy probability distribution.

from math import *
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

def read_E(filename):
    infile = open(filename, 'r')
    # Elements to be read in file:
    E = [];
    # Read lines except for the first one:
    lines = infile.readlines()
    first_line = lines[0]
    words_in_first_line = first_line.split()
    T = float(words_in_first_line[-1])
    for line in lines[2:-1]:
    	words = line.split()
        E.append(float(words[0]))
    final_line = lines[-1]
    E_variance = float(final_line.split()[-1])
    infile.close()
    return T, E, E_variance

# Fetching data by a call on read_E_and_M:
T, E, E_variance = read_E('Ising2D_d.txt')

print "Energy per spin std.: ", E_variance

# Plot results:
fig, ax = plt.subplots(1)
plt.rcParams.update({'font.size': 14})
n, bins, patches = plt.hist(E, 15, normed=1, facecolor='red', alpha=0.75, label='$\sigma_E^2 = %.3f$' % E_variance)
(mu, sigma) = norm.fit(E)
y = mlab.normpdf( bins, mu, sigma)
l = plt.plot(bins, y, 'k--', linewidth=2, label='Best fit Gaussian, $\sigma = %.2f$' % sigma)

# Checking the returned variables:
print "n: ", n
print "Bins: ", bins
print "Sum of bin values: ", np.sum(n)
print "Sum of bin areas: ", np.sum(n*np.diff(bins))

plt.xlabel('Energy per spin, $E/n_{\mathrm{spins}}$')
plt.ylabel('Probability, $P(E)$')
plt.legend(loc='upper right',fancybox='True')
plt.title('Energy histogram after reached equilibrium, T = %.1f in units of kT/J' % T)
plt.grid(True)
ax.set_ylim([0.0,4.0])
plt.savefig('Project4_d_probability_highT.eps', format='eps', dpi=1000)
plt.show()