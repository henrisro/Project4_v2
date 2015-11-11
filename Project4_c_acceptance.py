
# Makefile program for running Project4_c_acceptance.cpp

import os
os.system('g++ -c Project4_c_acceptance.cpp lib.cpp')
os.system('g++ -o Project4_c_acceptance.out Project4_c_acceptance.o lib.o -O3')
# Output file should be provided here as well:
os.system('./Project4_c_acceptance.out Ising2D_c_acceptance_ordered.txt')