// Project 4 e) 
// Plotting various quantities as functions of temperature to see indications
// of a phase transition. We study a periodic Ising model for sizes 20x20,
// 40x40, 60x60 and 80x80 in the temperature range T = [2.0,2.4]. 
// Exact transition occurs at T_c = 2.269 in the thermodynamic limit.

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "lib.h"
#include <time.h>
using namespace  std;

ofstream ofile;

// inline function for periodic boundary conditions
inline int periodic(int i, int limit, int add) { 
  return (i+limit+add) % (limit);
}
// Function to read in data from screen  
void read_input(int&, int&, double&, double&, double&);
// Function to initialise energy and magnetization
void initialize(int, double, int **, double&, double&);
// The metropolis algorithm 
void Metropolis(int, long&, int **, double&, double&, double *);
// prints to file the results of the calculations  
void output(int, int, double, double *, double);

int main(int argc, char* argv[])
{
  char *outfilename;
  long idum;
  int **spin_matrix, n_spins, mcs, mcs_i;
  double w[17], average[5], initial_temp, final_temp, E, M, temp_step;
  // This should really be adjusted as temperature is changed:
  int steady_state_tolerance_cycles = 5E3;
  double calculation_time;

  // Read in output file, abort if there are too few command-line arguments
  if( argc <= 1 ){
    cout << "Bad Usage: " << argv[0] << 
      " read also output file on same line" << endl;
    exit(1);
  }
  else{
    outfilename=argv[1];
  }
  ofile.open(outfilename);
  //    Read in initial values such as size of lattice, temp and cycles
  read_input(n_spins, mcs, initial_temp, final_temp, temp_step);
  int effective_mcs = mcs - steady_state_tolerance_cycles;
  int initial_temp_bol = 1; // Boolean variable: 1 if initial configuration.
  spin_matrix = (int**) matrix(n_spins, n_spins, sizeof(int));
  // Initialization of spin matrix (ordered initial state for low temps):
  E = M = 0.;
  double temperature = initial_temp;
  initialize(n_spins, temperature, spin_matrix, E, M);

  idum = -1; // random starting point
  for (double temperature = initial_temp; temperature <= final_temp; temperature+=temp_step){
    // setup array for possible energy changes
    for( int de =-8; de <= 8; de++) w[de+8] = 0;
    for( int de =-8; de <= 8; de+=4) w[de+8] = exp(-de/temperature);
    // initialise array for expectation values
    for( int i = 0; i < 5; i++) average[i] = 0.;

    // Manual initialization here not needed 
    // – use previous spin matrix as initial for the next computation.

    //initialize(n_spins, temperature, spin_matrix, E, M);
    // start Monte Carlo computation
    clock_t start, finish;
    start = clock();
    for (int cycles = 1; cycles <= mcs; cycles++){
      Metropolis(n_spins, idum, spin_matrix, E, M, w);
      // update expectation values
      // Initialize time:
      if (cycles >= steady_state_tolerance_cycles) {
        //cout << "Hit! Average contributions counted. cycles = " << cycles << endl;
        average[0] += E;    average[1] += E*E;
        average[2] += M;    average[3] += M*M; average[4] += fabs(M);
      }
    }
    finish = clock();
    calculation_time = (finish - start)/(double)CLOCKS_PER_SEC;
    // write final results to file:
    //cout << "Final effective number of cycles: " << effective_mcs << endl;
    output(n_spins, effective_mcs, temperature, average, calculation_time);
  }
  free_matrix((void **) spin_matrix); // free memory
  ofile.close();  // close output file
  return 0;
}

// read in input data
void read_input(int& n_spins, int& mcs, double& initial_temp, 
		double& final_temp, double& temp_step)
{
  cout << "Number of Monte Carlo trials: "; 
  cin >> mcs;
  cout << "Lattice size or number of spins (x and y equal): ";
  cin >> n_spins;
  cout << "Initial temperature with dimension energy: ";
  cin >> initial_temp;
  cout << "Final temperature with dimension energy: ";
  cin >> final_temp;
  cout << "Temperature step with dimension energy: ";
  cin >> temp_step;
} // end of function read_input


// function to initialise energy, spin matrix and magnetization
void initialize(int n_spins, double temperature, int **spin_matrix, 
		double& E, double& M)
{
  // setup spin matrix and intial magnetization
  for(int y =0; y < n_spins; y++) {
    for (int x= 0; x < n_spins; x++){
      spin_matrix[y][x] = 1; // spin orientation for the ground state
      M +=  (double) spin_matrix[y][x];
    }
  }

  // long idum_dum = -2;
  // for(int y =0; y < n_spins; y++) {
  //    for (int x= 0; x < n_spins; x++){
  //      double a = (double) ran1(&idum_dum);
  //      int r = 1;
  //      if (a < 0.5) {r = -1;}
  //      //cout << r;
  //      spin_matrix[y][x] = r; // spin orientation for the random state
  //      M += spin_matrix[y][x];
  //    }
  //    //cout << endl;
  // }
  // setup initial energy
  for(int y =0; y < n_spins; y++) {
    for (int x= 0; x < n_spins; x++){
      E -=  (double) spin_matrix[y][x]*
	(spin_matrix[periodic(y,n_spins,-1)][x] +
	 spin_matrix[y][periodic(x,n_spins,-1)]);
    }
  }
}// end function initialise

void Metropolis(int n_spins, long& idum, int **spin_matrix, double& E, double&M, double *w)
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
      }
    }
  }
} // end of Metropolis sampling over spins


void output(int n_spins, int mcs, double temperature, double *average, double calc_time)
{
  double norm = 1/((double) (mcs));  // divided by total number of cycles 
  double Eaverage = average[0]*norm;
  double E2average = average[1]*norm;
  double Maverage = average[2]*norm;
  double M2average = average[3]*norm;
  double Mabsaverage = average[4]*norm;
  // all expectation values are per spin, divide by 1/n_spins/n_spins
  double Evariance = (E2average - Eaverage*Eaverage)/n_spins/n_spins;
  double Mvariance = (M2average - Mabsaverage*Mabsaverage)/n_spins/n_spins;
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << setw(15) << setprecision(8) << temperature;
  ofile << setw(15) << setprecision(8) << Eaverage/n_spins/n_spins;
  ofile << setw(15) << setprecision(8) << Evariance/temperature/temperature;
  ofile << setw(15) << setprecision(8) << Maverage/n_spins/n_spins;
  ofile << setw(15) << setprecision(8) << Mvariance/temperature;
  ofile << setw(15) << setprecision(8) << Mabsaverage/n_spins/n_spins;
  ofile << setw(15) << setprecision(8) << calc_time << endl;;
} // end output function