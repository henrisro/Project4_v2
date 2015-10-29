//   Project 4 e) Create snapshot
//   Store spin matrix to visualize for equilibrium situation for different
//   temperatures.

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <time.h>
#include <random>
#include "lib.h"

using namespace  std;
ofstream ofile;

// inline function for periodic boundary conditions
inline int periodic(int i, int limit, int add) { 
  return (i+limit+add) % (limit);
}
// Function to read in data from screen  
void read_input(int&, double&, int&);
// Function to initialise energy and magnetization
void initialize(int, int **);
// The metropolis algorithm 
void Metropolis(int, long&, int **, double *);
// prints to file the results of the calculations  
void output(int, int **);

int main(int argc, char* argv[])
{
  char *outfilename;
  long idum;
  int **spin_matrix, n_spins, mcs;
  double w[17], average[5], temperature; //E, M;

  // Read in output file, abort if there are too few command-line arguments
  if( argc <= 1 ){
    cout << "Bad Usage: " << argv[0] << 
      " read also output file on same line" << endl;
    exit(1);
  }
  else{
    outfilename=argv[1];
  }
  read_input(n_spins, temperature, mcs);
  // Write header in output file:
  ofile.open(outfilename);
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << "Final number of Monte Carlo trials: " << mcs << ", and temperature: " << temperature << endl;
  ofile << " Spin matrix in equilibrium situation: " << endl;
  //    Read in initial values such as size of lattice, temp and cycles
  idum = -1; // random starting point
  spin_matrix = (int**) matrix(n_spins, n_spins, sizeof(int));
  // setup array for possible energy changes
  for( int de =-8; de <= 8; de++) w[de+8] = 0;
  for( int de =-8; de <= 8; de+=4) w[de+8] = exp(-de/temperature);
  // initialise array for expectation values
  //for( int i = 0; i < 5; i++) average[i] = 0.;
  initialize(n_spins, spin_matrix);
  // start Monte Carlo computation:
  for (int cycles = 1; cycles <= mcs; cycles++){
    Metropolis(n_spins, idum, spin_matrix, w);
    // update expectation values
  }
  // print results
  output(n_spins, spin_matrix);
  free_matrix((void **) spin_matrix); // free memory
  
  ofile.close();  // close output file
  return 0;
}

// read in input data
void read_input(int& n_spins, double& temperature, int& mcs)
{
  cout << "Final number of MC trials: ";
  cin >> mcs;
  cout << "Lattice size or number of spins (x and y equal): ";
  cin >> n_spins;
  cout << "Temperature with dimension energy: ";
  cin >> temperature;
} // end of function read_input


// function to initialise energy, spin matrix and magnetization
void initialize(int n_spins, int **spin_matrix)
{

  // Setup spin matrix and intial magnetization; all spins up.
  // Commment/uncomment this section depending on use
  // Ground state configuration gives good convergence for low temperatures

  for(int y =0; y < n_spins; y++) {
    for (int x= 0; x < n_spins; x++){
      spin_matrix[y][x] = 1; // spin orientation for the ground state
    }
  }
  
  // Setup spin matrix and intial magnetization; random start configuration.
  // Commment/uncomment this section depending on use.
  // Random spin configuration gives good convergence for high temperatures

  // long idum_dum = 2;
  // for(int y =0; y < n_spins; y++) {
  //    for (int x= 0; x < n_spins; x++){
  //      double a = (double) ran1(&idum_dum);
  //      int r = 1;
  //      if (a < 0.5) {r = -1;}
  //      spin_matrix[y][x] = r; // spin orientation for the random state
  //      M +=  (double) spin_matrix[y][x];
  //    }
  //  }

}// end function initialise

void Metropolis(int n_spins, long& idum, int **spin_matrix, double *w)
{
  // loop over all spins
  for(int y =0; y < n_spins; y++) {
    for (int x= 0; x < n_spins; x++){
      int ix = (int) (ran1(&idum)*(double)n_spins);
      int iy = (int) (ran1(&idum)*(double)n_spins);
      int deltaE =  2*spin_matrix[iy][ix]*
  (spin_matrix[iy][periodic(ix,n_spins,-1)]+
   spin_matrix[periodic(iy,n_spins,-1)][ix] +
   spin_matrix[iy][periodic(ix,n_spins,1)] +
   spin_matrix[periodic(iy,n_spins,1)][ix]);
      if ( ran1(&idum) <= w[deltaE+8] ) {
  spin_matrix[iy][ix] *= -1;  // flip one spin and accept new spin config
      }
    }
  }
} // end of Metropolis sampling over spins

void output(int n_spins, int **spin_matrix)
{
  // Write final spin matrix to file:
  for (int nx = 0; nx < n_spins; nx++) {
      for (int ny = 0; ny < n_spins; ny++) {
          ofile << setw(5) << setprecision(2) << spin_matrix[ny][nx];
      }
      ofile << endl;
  }
} // end output function