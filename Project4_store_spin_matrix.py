
# Makefile program for running Project4_store_spin_matrix.cpp and
# plot a snapshot of the equilibrium state. This is done by
# visulizing the binary values of the spin matrix.

import os
os.system('g++ -c Project4_store_spin_matrix.cpp lib.cpp')
os.system('g++ -o Project4_store_spin_matrix.out Project4_store_spin_matrix.o lib.o -O2')
# Output file should be provided here as well:
os.system('./Project4_store_spin_matrix.out Spin_matrix.txt')

from math import *
import numpy as np
import matplotlib.pyplot as plt

def spin_matrix_plot(filename):
    infile = open(filename, 'r')

    plt.rcParams.update({'font.size': 14})
    fig, ax = plt.subplots(1)
    # Read lines except for the first one:
    lines = infile.readlines()
    second_line = lines[2]
    first_line = lines[0]
    n_spins = len(second_line.split())
    T = float(first_line.split()[-1])
    ny = 0
    for line in lines[2:]:
    	words = line.split()
    	y_coor = ny/float(n_spins)
    	for nx in range(n_spins):
    		x_coor = nx/float(n_spins)
    		# Plot spin up and down in different colors:
    		if int(words[nx]) == +1:
    			ax.plot([x_coor], [y_coor], 'wo')
    		else: ax.plot([x_coor], [y_coor], 'ko')
        ny += 1
    ax.set_xlabel('$x$ coordinate in grid')
    ax.set_ylabel('$y$ coordinate in grid')
    ax.set_title('Snapshot of equilibirum state in %dx%d grid. T = %.1f in units of kT/J' % (n_spins,n_spins,T)) 
    ax.set_xlim([-0.02,1.0])
    ax.set_ylim([-0.02,1.0])
    infile.close()
    plt.savefig('Snapshot_T25.eps', format='eps', dpi=1000)
    plt.show()
    return 0

# Visualize by calling function:
spin_matrix_plot('Spin_matrix.txt')