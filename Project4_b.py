# Makefile program for running Project4_b.cpp

import os

os.system('g++ -c Project4_b.cpp lib.cpp')
os.system('g++ -o Project4_b.out Project4_b.o lib.o -O2')
# Output file should be provided here as well:
os.system('./Project4_b.out Ising2D_b.txt')