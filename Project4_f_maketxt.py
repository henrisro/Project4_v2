
# Makefile program for running Project4_e.cpp

import os
os.system('g++ -c Project4_e.cpp lib.cpp')
os.system('g++ -o Project4_e.out Project4_e.o lib.o -O3')
# Output file should be provided here as well:
os.system('./Project4_e.out Ising2D_f_L20_time_comparison.txt')