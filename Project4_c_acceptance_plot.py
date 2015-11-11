# Read data from 'Ising2D_c_acceptance_<ordered/disordered>.txt' and plot energy and abs(magnetization) as
# functions of the number of MC cycles.

from math import *
import numpy as np
import matplotlib.pyplot as plt

def read_acceptance(filename):
    infile = open(filename, 'r')
    # Elements to be read in file:
    T = []; E = []; M = []; accepted = []
    # Read lines except for the first one:
    lines = infile.readlines()
    first_line = lines[0]
    words_in_first_line = first_line.split()
    N = float(words_in_first_line[-1])
    for line in lines[2:]:
    	words = line.split()
        T.append(float(words[0]))
        E.append(float(words[1]))
        M.append(float(words[2]))
        accepted.append(log10(float(words[3])))
    infile.close()
    return T, E, M, accepted, N

# Number of Monte Carlo cycles
# Fetching data by a call on read_E_and_M:
T1, E1, M1, accepted1, N1 = read_acceptance('Ising2D_c_acceptance_ordered.txt')
T2, E2, M2, accepted2, N2 = read_acceptance('Ising2D_c_acceptance_disordered.txt')

plt.rcParams.update({'font.size': 14})
fig, ax = plt.subplots(1)
ax.plot(T1,accepted1,label='Ordered, $N_{MC}=%g$' % N1)
ax.plot(T2,accepted2,label='Disordered, $N_{MC}=%g$' % N2)
ax.set_xlabel('Temperature, $T$')
ax.set_ylabel('Log of total fraction of accepted configurations, $\log_{10} A(N_{MC})$')
ax.set_title('Total fraction of accepted configurations') 
ax.grid()
ax.legend(loc='best',fancybox='True',shadow='True')
#ax.set_ylim([-2.00,0.0])
plt.savefig('Project4_c_acceptance.eps', format='eps', dpi=1000)
plt.show()