//   Project 4 c)
//   Program to study various expectation values as functions
//   of the number of Monte Carlo cycles. We use a 20x20 grid Ising model
//   with periodic boundary conditions. Values are plotted with python script. 

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <time.h>
#include <random>
#include "lib.h"

using namespace std;
ofstream ofile;

// inline function for periodic boundary conditions
inline int periodic(int i, int limit, int add) { 
  return (i+limit+add) % (limit);
}
// Function to read in data from screen  
void read_input(int&, double&, int&);
// Function to initialise energy and magnetization
void initialize(int, double, int **, double&, double&);
// The metropolis algorithm 
void Metropolis(int, long&, int **, double&, double&, double *, int&);
// prints to file the results of the calculations  
void output(int, int, double *, int);

int main(int argc, char* argv[])
{
  char *outfilename;
  long idum;
  int **spin_matrix, n_spins, mcs, mcs_i;
  double w[17], average[5], temperature, E, M;
  int no_accepted = 0;

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
  ofile << "     Energy:      abs(Magnetization):    No of accepted moves: " << endl;
  // Read in initial values such as size of lattice, temp and cycles
  idum = -3; // random starting point
  spin_matrix = (int**) matrix(n_spins, n_spins, sizeof(int));
  //    initialise energy and magnetization 
  E = M = 0.;
  // setup array for possible energy changes
  for( int de =-8; de <= 8; de++) w[de+8] = 0;
  for( int de =-8; de <= 8; de+=4) w[de+8] = exp(-de/temperature);
  // initialise array for expectation values
  for( int i = 0; i < 5; i++) average[i] = 0.;
  initialize(n_spins, temperature, spin_matrix, E, M);
  // Printing initial values of E and M to the result file:
  double initial_Eaverage = E/n_spins/n_spins;
  double initial_Mabsaverage = fabs(M)/n_spins/n_spins;
  ofile << setw(15) << setprecision(8) << initial_Eaverage;
  ofile << setw(15) << setprecision(8) << initial_Mabsaverage; 
  ofile << setw(15) << setprecision(8) << 0 << endl;
  // start Monte Carlo computation:
  for (int cycles = 1; cycles <= mcs; cycles++){
    mcs_i = cycles; // Current number of cycles.
    Metropolis(n_spins, idum, spin_matrix, E, M, w, no_accepted);
    // update expectation values
    average[0] += E;    average[1] += E*E;
    average[2] += M;    average[3] += M*M;   average[4] += fabs(M);
    output(n_spins, mcs_i, average, no_accepted);
  }
  // print results
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
void initialize(int n_spins, double temperature, int **spin_matrix, 
                double& E, double& M)
{

  // Setup spin matrix and intial magnetization; all spins up.
  // Commment/uncomment this section depending on use
  // Ground state configuration gives good convergence for low temperatures

  // for(int y =0; y < n_spins; y++) {
  //   for (int x= 0; x < n_spins; x++){
  //     spin_matrix[y][x] = 1; // spin orientation for the ground state
  //     M +=  (double) spin_matrix[y][x];
  //   }
  // }
  
  // Setup spin matrix and intial magnetization; random start configuration.
  // Commment/uncomment this section depending on use.
  // Random spin configuration gives good convergence for high temperatures

  long idum_dum = -2;
  for(int y =0; y < n_spins; y++) {
     for (int x= 0; x < n_spins; x++){
       double a = (double) ran1(&idum_dum);
       int r = 1;
       if (a < 0.5) {r = -1;}
       //cout << r;
       spin_matrix[y][x] = r; // spin orientation for the random state
       M += spin_matrix[y][x];
     }
     //cout << endl;
  }
  cout << "Initial magnetization: " << M << endl;

  // setup initial energy:
  for(int y =0; y < n_spins; y++) {
    for (int x= 0; x < n_spins; x++){
      E -=  (double) spin_matrix[y][x]*
  (spin_matrix[periodic(y,n_spins,-1)][x] +
   spin_matrix[y][periodic(x,n_spins,-1)]);
    }
  }
}// end function initialise

void Metropolis(int n_spins, long& idum, int **spin_matrix, double &E, double &M, double *w, int &acceptance)
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
        M += (double) 2*spin_matrix[iy][ix];
        E += (double) deltaE;
        acceptance += 1;
      }
    }
  }
} // end of Metropolis sampling over spins

void output(int n_spins, int mcs_i, double *average, int accepted)
{
  double norm = 1/((double) (mcs_i));  // divide by total number of cycles 
  double norm2 = 1.0/(n_spins*n_spins);   // divide by the total number of spins
  double Eaverage = average[0]*norm;
  double Mabsaverage = average[4]*norm;
  // all expectation values are per spin, divide by 1/n_spins/n_spins
  ofile << setw(15) << setprecision(8) << Eaverage*norm2;
  ofile << setw(15) << setprecision(8) << Mabsaverage*norm2;
  ofile << setw(15) << setprecision(8) << accepted << endl;
} // end output function