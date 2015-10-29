
# Makefile program for running Project4_c.cpp

import os
os.system('g++ -c Project4_c.cpp lib.cpp')
os.system('g++ -o Project4_c.out Project4_c.o lib.o -O2')
# Output file should be provided here as well:
os.system('./Project4_c.out Ising2D_c_disordered_highT.txt')