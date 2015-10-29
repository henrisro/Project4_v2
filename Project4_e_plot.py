# Read data from all the different txt-files and plot different quantities
# as function of temperature

from math import *
import numpy as np
import matplotlib.pyplot as plt

def read_quantities(filename):
    infile = open(filename, 'r')
    # Elements to be read in file:
    T = []; E = []; Cv = []; Chi = []; absM = [] 
    # Read lines except for the first one:
    lines = infile.readlines()
    for line in lines:
    	words = line.split()
        T.append(float(words[0]))
        E.append(float(words[1]))
        Cv.append(float(words[2]))
        # Skip the storage of M here.
        Chi.append(float(words[4]))
        absM.append(float(words[5]))
    infile.close()
    return T, E, Cv, Chi, absM

# Some exact results provided by Onsager (1944):
Tc_exact = 2.0/log(1.0+sqrt(2.0))

def magnetization(T):
    M = 0.0
    if T < Tc_exact:
       M = (1.0-(sinh(2/T))**(-4))**(1.0/8.0)
    return M

# Complete elliptic integrals of the first and second kind:
def K_minus(theta,k1):
    return 1.0/sqrt(1.0-k1**2*sin(theta)**2)

def K_plus(theta,k1):
    return sqrt(1.0-k1**2*sin(theta)**2)

def integral_K_minus(k1,N=401):
    # Evaluating complete elliptic integral numerically:
    theta = np.linspace(0,pi/2.0,N)
    h = theta[1]-theta[0]
    integral = 0.0
    # Simple integration with the trapezoidal rule:
    for i in range(N-1):
        integral += h/2.0*(K_minus(theta[i],k1)+K_minus(theta[i+1],k1))
    return integral

def integral_K_plus(k1,N=401):
    # Evaluating complete elliptic integral numerically:
    theta = np.linspace(0,pi/2.0,N)
    h = theta[1]-theta[0]
    integral = 0.0
    # Simple integration with the trapezoidal rule:
    for i in range(N-1):
        integral += h/2.0*(K_plus(theta[i],k1)+K_plus(theta[i+1],k1))
    return integral

def heat_capacity_per_spin(T):
    k1 = 2.0*sinh(2.0/T)/(cosh(2.0/T)*cosh(2.0/T))
    k1_doubleprime = 2*(tanh(2.0/T))**2-1.0
    k1_tripleprime = 1.0-(tanh(2.0/T))**2
    term1 = integral_K_minus(k1,N=401)
    term2 = integral_K_plus(k1,N=401)
    coth = cosh(2.0/T)/sinh(2.0/T)
    return 4.0/pi*(1.0/T*coth)**2*( term1 - term2 - k1_tripleprime*(pi/2.0+k1_doubleprime*term1) )

def internal_energy_per_spin(T):
    k1 = 2.0*sinh(2.0/T)/(cosh(2.0/T)*cosh(2.0/T))
    k1_doubleprime = 2*(tanh(2.0/T))**2-1.0
    term1 = integral_K_minus(k1,N=401)
    C = -cosh(2.0/T)/sinh(2.0/T)
    return C*( 1.0 + 2.0/pi*k1_doubleprime*term1 )

n = 201
T_vals = np.linspace(2.0,2.4,n)
M_exact = np.zeros(n)
E_exact = np.zeros(n)
Cv_exact = np.zeros(n)
for i in range(n):
    T_i = T_vals[i]
    M_exact[i] = magnetization(T_i)
    E_exact[i] = internal_energy_per_spin(T_i)
    Cv_exact[i] = heat_capacity_per_spin(T_i)

# Fetching data by a call on read_quantities:
T2, E2, Cv2, Chi2, absM2 = read_quantities('Ising2D_e_L20.txt')
T4, E4, Cv4, Chi4, absM4 = read_quantities('Ising2D_e_L40.txt')
T6, E6, Cv6, Chi6, absM6 = read_quantities('Ising2D_e_L60.txt')
T8, E8, Cv8, Chi8, absM8 = read_quantities('Ising2D_e_L80.txt')

plt.rcParams.update({'font.size': 14})
fig, ax = plt.subplots(1)
ax.plot(T2, E2,label='$L=20$')
ax.plot(T4, E4,label='$L=40$')
ax.plot(T6, E6,label='$L=60$')
ax.plot(T8, E8,label='$L=80$')
ax.plot([Tc_exact,Tc_exact],[min(E8),max(E8)],'m--',label='$T_c$ $\mathrm{exact}$')
ax.plot(T_vals, E_exact, 'k-', label='$\mathrm{Onsager}$ $\mathrm{solution}$')
ax.legend(loc='best',fancybox='True',shadow='True')
ax.set_xlabel('Temperature, $T$, in units of $kT/J$')
ax.set_ylabel('Average energy per spin, $\langle E \\rangle / n_{\mathrm{spins}}$')
ax.set_title('Average energy as function of temperature') 
ax.grid()
#plt.savefig('Project4_e_energy.eps', format='eps', dpi=1000)
plt.show()

fig, ax = plt.subplots(1)
ax.plot(T2, Cv2,label='$L=20$')
ax.plot(T4, Cv4,label='$L=40$')
ax.plot(T6, Cv6,label='$L=60$')
ax.plot(T8, Cv8,label='$L=80$')
ax.plot([Tc_exact,Tc_exact],[min(Cv_exact),max(Cv_exact)],'m--',label='$T_c$ $\mathrm{exact}$')
ax.plot(T_vals, Cv_exact, 'k-', label='$\mathrm{Onsager}$ $\mathrm{solution}$')
ax.set_xlim([2.0,2.4])
ax.set_ylim([0.6,3.0])
ax.legend(loc='best',fancybox='True',shadow='True')
ax.set_xlabel('Temperature, $T$, in units of $kT/J$')
ax.set_ylabel('Heat capacity per spin, $C_v / n_{\mathrm{spins}}$')
ax.set_title('Heat capacity as function of temperature') 
ax.grid()
#plt.savefig('Project4_e_heatcapacity.eps', format='eps', dpi=1000)
plt.show()

fig, ax = plt.subplots(1)
ax.plot(T2, absM2,label='$L=20$')
ax.plot(T4, absM4,label='$L=40$')
ax.plot(T6, absM6,label='$L=60$')
ax.plot(T8, absM8,label='$L=80$')
ax.plot([Tc_exact,Tc_exact],[min(absM8),max(absM8)],'m--',label='$T_c$ $\mathrm{exact}$')
ax.plot(T_vals, M_exact, 'k-', label='$\mathrm{Onsager}$ $\mathrm{solution}$')
ax.legend(loc='best',fancybox='True',shadow='True')
ax.set_xlabel('Temperature, $T$, in units of $kT/J$')
ax.set_ylabel('Average magnetization per spin, $\langle \mid \mathcal{M} \mid \\rangle / n_{\mathrm{spins}}$')
ax.set_title('Average magnetization as function of temperature') 
ax.grid()
#plt.savefig('Project4_e_magnetization.eps', format='eps', dpi=1000)
plt.show()

fig, ax = plt.subplots(1)
ax.plot(T2, Chi2,label='$L=20$')
ax.plot(T4, Chi4,label='$L=40$')
ax.plot(T6, Chi6,label='$L=60$')
ax.plot(T8, Chi8,label='$L=80$')
ax.plot([Tc_exact,Tc_exact],[min(Chi8),max(Chi8)],'m--',label='$T_c$ $\mathrm{exact}$')
ax.legend(loc='best',fancybox='True',shadow='True')
ax.set_xlabel('Temperature, $T$, in units of $kT/J$')
ax.set_ylabel('Susceptibility per spin, $\chi / n_{\mathrm{spins}}$')
ax.set_title('Susceptibility as function of temperature') 
ax.grid()
#plt.savefig('Project4_e_susceptibility.eps', format='eps', dpi=1000)
plt.show()