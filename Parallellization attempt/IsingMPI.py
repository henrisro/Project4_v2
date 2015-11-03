# Makefile program for running isingMPI.cpp

import os

os.system('mpic++ -o ising.x IsingMPI.cpp -O3')
# Output file should be provided here as well:
os.system('./ising.x Ising2D_f_L20_time_comparison_MPI.txt')