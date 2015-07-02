#include <iostream>
#include "nbody.h"
#include "nbodyio.h"

using namespace std;

/*-----------------------------------------------------------------------------
 *  get_snapshot  --  reads a single snapshot from the input stream cin.
 *                    Only the particle data is read in- the main program is
 *                    responsible for reading in particle number and time.
 *                    The system is normalized with the center of mass at the
 *                    origin and zero net momentum.
 *-----------------------------------------------------------------------------
 */

void get_snapshot(real mass[], real pos[][NDIM], real vel[][NDIM], int n){

	real (*cmass) = new real[NDIM];      //center of mass
	real (*moment) = new real[NDIM];     //net momentum
	for(int k = 0; k < NDIM; k++){
		cmass[k] = 0;
		moment[k] = 0;
	}

	real tmass = 0;
	for(int i = 0; i < n; i++){
		real m;
		cin >> m;
		mass[i] = G*m;                  // mass of particle i
		tmass += m;
		for(int k = 0; k < NDIM; k++){
			real p;
			cin >> p;
			pos[i][k] = p;              // position of particle i
			cmass[k] += p*m;
		}
		for(int k = 0; k < NDIM; k++){
			real v;
			cin >> v;
			vel[i][k] = v;              // velocity of particle i
			moment[k] += v*m;
		}
	}

	for(int k = 0; k < NDIM; k++){
		real offset = cmass[k]/tmass;
		real velocity = moment[k]/tmass;
		for(int i = 0; i < n; i++){
			pos[i][k] -= offset;
			vel[i][k] -= velocity;
		}
	}
}

/*-----------------------------------------------------------------------------
 *  put_snapshot  --  writes a single snapshot on the output stream cout.
 *  note: unlike get_snapshot(), put_snapshot handles particle number and time
 *-----------------------------------------------------------------------------
 */

void put_snapshot(const real mass[], const real pos[][NDIM],
				  const real vel[][NDIM], const real dst[],
				  int n, real t){

	cout.precision(16);
	cout << n << ' ' << t << endl;
	for(int i = 0; i < n; i++){
		cout << (mass[i] / G);
		for(int k = 0; k < NDIM; k++){ cout << ' ' << pos[i][k]; }
		for(int k = 0; k < NDIM; k++){ cout << ' ' << vel[i][k]; }
		cout << endl;
	}

	int d = n*(n-1)/2;
	for(int i = 0; i < d; i++){ cout << dst[i] << ' '; }
	cout << endl;
}

/*-----------------------------------------------------------------------------
 *  write_diagnostics  --  writes diagnostics on the error stream cerr:
 *                         current time; number of steps so far; KE, PE, and
 *                         total energy; absolute and relative energy errors
 *                         since the start of the run.
 *                         If x_flag (for eXtra data) is true, all internal
 *                         data is dumped for each particle (mass, position,
 *                         velocity, acceleration, and jerk).
 *
 *  note: the kinetic energy is calculated here, while the potential energy is
 *        calculated in the function get_acc_jrk_pot_coll().
 *-----------------------------------------------------------------------------
 */

void write_diagnostics(const real mass[], const real pos[][NDIM],
					   const real vel[][NDIM], const real acc[][NDIM],
					   const real jrk[][NDIM], int n, real t, real epot,
					   int nsteps, bool x_flag){

	static real einit = 0;

	real ekin = 0;                       // kinetic energy of the n-body system
	for(int i = 0; i < n; i++){
		for(int k = 0; k < NDIM; k++){
			ekin += 0.5 * mass[i] * vel[i][k] * vel[i][k];
		}
	}

	real etot = (ekin + epot)/G;         // total energy of the n-body system

	if(einit == 0){                      // at first pass, record the initial energy
		einit = etot;
	}

	cerr << "at time t = " << t << " , after " << nsteps
		 << " steps :\n  E_kin = " << ekin
		 << " , E_pot = " << epot
		 << " , E_tot = " << etot << endl;
	cerr << "                "
		 << "absolute energy error: E_tot - E_init = "
		 << etot - einit << endl;
	cerr << "                "
		 << "relative energy error: (E_tot - E_init) / E_init = "
		 << (etot - einit) / einit << endl;

	if(!x_flag){ return; }
	cerr << "  Internal data: \n";
	for(int i = 0; i < n; i++){
		cerr << "    Data for particle " << i+1 << " : " << endl;
		cerr << "      Mass: ";
		cerr << mass[i];
		cerr << "\n      Pos:  ";
		for(int k = 0; k < NDIM; k++)
			cerr << ' ' << pos[i][k];
		cerr << "\n      Vel:  ";
		for(int k = 0; k < NDIM; k++)
			cerr << ' ' << vel[i][k];
		cerr << "\n      Acc:  ";
		for(int k = 0; k < NDIM; k++)
			cerr << ' ' << acc[i][k];
		cerr << "\n      jrk: ";
		for(int k = 0; k < NDIM; k++)
			cerr << ' ' << jrk[i][k];
		cerr << endl;
	}
}