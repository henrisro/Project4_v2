# Read data from 'Ising2D_c_<ordered/disordered>.txt' and plot energy and abs(magnetization) as
# functions of the number of MC cycles.

from math import *
import numpy as np
import matplotlib.pyplot as plt

def read_E_and_M(filename):
    infile = open(filename, 'r')
    # Elements to be read in file:
    E = []; M = []; accepted = []
    # Read lines except for the first one:
    lines = infile.readlines()
    first_line = lines[0]
    words_in_first_line = first_line.split()
    T = float(words_in_first_line[-1])
    for line in lines[2:]:
    	words = line.split()
        E.append(float(words[0]))
        M.append(float(words[1]))
        accepted.append(int(words[2]))
    infile.close()
    return T, E, M, accepted

# Fetching data by a call on read_E_and_M:
T1, E1, M1, accepted1 = read_E_and_M('Ising2D_c_ordered_highT.txt')
T2, E2, M2, accepted2 = read_E_and_M('Ising2D_c_disordered_highT.txt')

plt.rcParams.update({'font.size': 13})
fig, ax = plt.subplots(1)
ax.plot(E1,label='Ordered initial state')
ax.plot(E2,label='Disordered initial state')
ax.set_xlabel('Number of Monte Carlo cycles, $N_{MC}$')
ax.set_ylabel('Average energy per spin, $\langle E \\rangle / n_{\mathrm{spins}}$')
ax.set_title('Average energy as function of number of MC cycles. \n T = %.1f in units of kT/J' % T1) 
ax.grid()
ax.legend(loc='best',fancybox='True',shadow='True')
ax.set_ylim([-2.00,0.0])
plt.savefig('Project4_c_energy_highT.eps', format='eps', dpi=1000)
plt.show()

fig2, ax2 = plt.subplots(1)
ax2.plot(M1,label='Ordered initial state')
ax2.plot(M2,label='Disordered initial state')
ax2.set_xlabel('Number of Monte Carlo cycles, $N_{MC}$')
ax2.set_ylabel('Average magnetization per spin, $\langle \mid \mathcal{M} \mid \\rangle  / n_{\mathrm{spins}}$')
ax2.set_title('Average magnetization as function of number of MC cycles. \n T = %.1f in units of kT/J' % T1) 
ax2.grid()
ax2.legend(loc='best',fancybox='True',shadow='True')
ax2.set_ylim([0.0,1.05])
plt.savefig('Project4_c_magnetization_highT.eps', format='eps', dpi=1000)
plt.show()

# fig3, ax3 = plt.subplots(1)
# ax3.plot(accepted,'g-')
# ax3.set_xlabel('Number of Monte Carlo cycles, $N_{MC}$')
# ax3.set_ylabel('Number of accepted configurations')
# ax3.set_title('Accepted configurations as function of MC cycles. \n T = %.1f in units of kT/J' % T)
# ax3.grid()
# #plt.savefig('Project4_c_accepted_ordered.eps', format='eps', dpi=1000)
# plt.show()