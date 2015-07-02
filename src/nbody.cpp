/*=============================================================================
 *
 *  nbody.cpp: an N-body integrator with a variable global time step,
 *             using the Hermite integration scheme.
 *
 *           ref.: Hut, P., Makino, J. & McMillan, S., 1995,
 *                  Astrophysical Journal Letters 443, L93-L96.
 *_____________________________________________________________________________
 *
 *  External data format:
 *
 *     The program expects input of a single snapshot of an N-body system,
 *     in the following format: the number of particles in the snapshot n;
 *     the time t; mass mi, position ri and velocity vi for each particle i,
 *     with position and velocity given through their three Cartesian
 *     coordinates, divided over separate lines as follows:
 *
 *                      n t
 *                      m1 r1_x r1_... v1_x v1_...
 *                      m2 r2_x r2_... v2_x v2_...
 *                      ...
 *                      mn rn_x rn_... vn_x vn_...
 *
 *     Output is written in the same format, so that old output can be used to
 *     restart a simulation where it left off, with the addition of one more
 *     line that's ignored on input. This is a list of all distances between
 *     pairs of particles, included for convenience, so that graphing software
 *     need not spend time recalculating it.
 *
 *  Internal data format:
 *
 *     The data for an N-body system is stored internally as a 1-dimensional
 *     array for the masses, and 2-dimensional arrays for the positions,
 *     velocities, accelerations and jrks of all particles.
 */

#include <iostream>
#include <cstdlib>    // for atoi() and atof()
#include <unistd.h>   // for getopt()
#include "nbody.h"
#include "nbodyio.h"
#include "evolve.h"

using namespace std;

void evolve(const real mass[], real pos[][NDIM], real vel[][NDIM], real dst[],
			int n, real t, real dt_param, real dt_dia, real dt_out, real dt_tot,
			bool x_flag);

bool read_options(int argc, char *argv[], real & dt_param, real & dt_dia,
				  real & dt_out, real & dt_tot, bool & x_flag);

/*-----------------------------------------------------------------------------
 *  main  --  reads options, reads a snapshot, and launches the integrator
 *-----------------------------------------------------------------------------
 */

int main(int argc, char *argv[]){
	real  dt_param = 0.03;     // control parameter to determine time step size
	real  dt_dia = 0;          // time interval between diagnostics output
	real  dt_out = 60;         // time interval between output of snapshots
	real  dt_tot = 3600;       // duration of the integration; default 1 hour
	bool  x_flag = false;      // if true: extra debugging diagnostics output

	if(!read_options(argc, argv, dt_param, dt_dia, dt_out, dt_tot, x_flag)){
		return 1;                // halt criterion detected by read_options()
	}

	int n;                       // number of particles in the N-body system
	cin >> n;

	real t;                      // time
	cin >> t;

	real *mass = new real[n];                  // masses for all particles
	real (*pos)[NDIM] = new real[n][NDIM];     // positions for all particles
	real (*vel)[NDIM] = new real[n][NDIM];     // velocities for all particles
	real (*dst) = new real[n*(n-1)/2];         // distances between all particle pairs

	get_snapshot(mass, pos, vel, n);

	cerr << "Starting a Hermite integration for a " << n
		 << "-body system,\n  from time t = " << t
		 << " with time step control parameter dt_param = " << dt_param
		 << "  until time " << t + dt_tot
		 << " ,\n  with diagnostics output interval dt_dia = "
		 << dt_dia << ",\n  and snapshot output interval dt_out = "
		 << dt_out << "." << endl;

	evolve(mass, pos, vel, dst, n, t,
		   dt_param, dt_dia, dt_out, dt_tot, x_flag);

	delete[] mass;
	delete[] pos;
	delete[] vel;
	delete[] dst;
}

/*-----------------------------------------------------------------------------
 *  read_options  --  reads the command line options.
 *  note: when the help option -h is invoked, or an unknown option encountered,
 *        the return value is set to false to prevent further execution.
 *-----------------------------------------------------------------------------
 */

bool read_options(int argc, char *argv[], real & dt_param, real & dt_dia,
				  real & dt_out, real & dt_tot, bool & x_flag){
	int c;
	while((c = getopt(argc, argv, "ha:e:o:t:x")) != -1){
		switch(c){
			case 'a': dt_param = atof(optarg);
					  break;
			case 'd': dt_dia = atof(optarg);
					  break;
			case 'o': dt_out = atof(optarg);
					  break;
			case 't': dt_tot = atof(optarg);
					  break;
			case 'x': x_flag = true;
					  break;
			case 'h': // fallthrough
			case '?': cerr << "usage: " << argv[0]
						   << " [-h (for help)]"
						   << " [-a step size control parameter]\n"
						   << "         [-d diagnostics interval]"
						   << " [-o output interval]\n"
						   << "         [-t total duration]"
						   << " [-x (extra debugging diagnostics)]"
						   << endl;
					  return false; // execution should stop after help or error
			}
	}

	return true; // continue program execution
}

/*-----------------------------------------------------------------------------
 *  evolve  --  integrates an N-body system, for a total duration dt_tot.
 *              Snapshots are sent to the standard output stream once every
 *              time interval dt_out.  Diagnostics are sent to the standard
 *              error stream once every time interval dt_dia.
 *
 *  In order to write diagnostics, we first have to calculate the potential
 *  energy with get_acc_jrk_pot_coll(). This also gives initial accelerations,
 *  jerks, and the collision time scale estimate, all of which are needed
 *  before we can enter the main integration loop.
 *
 *  note: the integration time step is global, but variable. Before each step
 *        we use the collision time estimate multiplied by dt_param (the
 *        accuracy parameter) to obtain the new time step size.
 *-----------------------------------------------------------------------------
 */

void evolve(const real mass[], real pos[][NDIM], real vel[][NDIM], real dst[],
			int n, real t, real dt_param, real dt_dia, real dt_out, real dt_tot,
			bool x_flag){

	real (* acc)[NDIM] = new real[n][NDIM];  // accelerations and jerks
	real (* jrk)[NDIM] = new real[n][NDIM];  // for all particles

	real (* old_pos)[NDIM] = new real[n][NDIM];
	real (* old_vel)[NDIM] = new real[n][NDIM];
	real (* old_acc)[NDIM] = new real[n][NDIM];
	real (* old_jrk)[NDIM] = new real[n][NDIM];

	real epot;                // potential energy of the n-body system
	real coll_time;           // collision (close encounter) time scale

	get_acc_jrk_pot_coll(mass, pos, vel, acc, jrk, dst, n, epot, coll_time);

	write_diagnostics(mass, pos, vel, acc, jrk,
					  n, t, epot, 0, x_flag);

	put_snapshot(mass, pos, vel, dst, n, t);

	real t_dia = t + dt_dia;  // next time for diagnostics output
	real t_out = t + dt_out;  // next time for snapshot output
	real t_end = t + dt_tot;  // final time, to finish the integration

	int nsteps = 0;           // number of integration time steps completed
	while(t < t_end){
		real dt = dt_param * coll_time;
		evolve_step(mass, pos, vel, acc, jrk, dst,
					old_pos, old_vel, old_acc, old_jrk,
					n, dt, epot, coll_time);
		t += dt;
		nsteps++;
		if(dt_dia > 0 && t >= t_dia){
			write_diagnostics(mass, pos, vel, acc, jrk,
							  n, t, epot, nsteps, x_flag);
			do{ t_dia += dt_dia; } while(t_dia < t);
		}
		if(t >= t_out){
			put_snapshot(mass, pos, vel, dst, n, t);
			do{ t_out += dt_out; } while(t_out < t);
		}
	}

	if(dt_dia == 0 || t > (t_dia - dt_dia)){
		write_diagnostics(mass, pos, vel, acc, jrk,
						  n, t, epot, nsteps, x_flag);
	}

	delete[] acc;
	delete[] jrk;
	delete[] old_pos;
	delete[] old_vel;
	delete[] old_acc;
	delete[] old_jrk;
}