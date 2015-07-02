#include <cmath>      // to include sqrt(), etc.
#include <cfloat>     // for DBL_MAX
#include "nbody.h"
#include "evolve.h"

using namespace std;

/*-----------------------------------------------------------------------------
 *  evolve_step  --  takes one integration step for an N-body system, using the
 *                   Hermite algorithm.
 *-----------------------------------------------------------------------------
 */

void evolve_step(const real mass[],
				 real pos[][NDIM], real vel[][NDIM],
				 real acc[][NDIM], real jrk[][NDIM], real dst[],
				 real old_pos[][NDIM], real old_vel[][NDIM],
				 real old_acc[][NDIM], real old_jrk[][NDIM],
				 int n, real dt, real & epot, real & coll_time){

	for(int i = 0; i < n; i++){
		for(int k = 0; k < NDIM; k++){
		  old_pos[i][k] = pos[i][k];
		  old_vel[i][k] = vel[i][k];
		  old_acc[i][k] = acc[i][k];
		  old_jrk[i][k] = jrk[i][k];
		}
	}

	predict_step(pos, vel, acc, jrk, n, dt);
	get_acc_jrk_pot_coll(mass, pos, vel, acc, jrk, dst, n, epot, coll_time);
	correct_step(pos, vel, acc, jrk, old_pos, old_vel, old_acc, old_jrk, n, dt);
}

/*-----------------------------------------------------------------------------
 *  predict_step  --  takes the first approximation of one Hermite integration
 *                    step, advancing the positions and velocities through a
 *                    Taylor series development up to the order of the jrks.
 *-----------------------------------------------------------------------------
 */

void predict_step(real pos[][NDIM], real vel[][NDIM],
				  const real acc[][NDIM], const real jrk[][NDIM],
				  int n, real dt){
	for(int i = 0; i < n; i++){
		for(int k = 0; k < NDIM; k++){
			pos[i][k] += vel[i][k]*dt + acc[i][k]*dt*dt/2
									  + jrk[i][k]*dt*dt*dt/6;
			vel[i][k] += acc[i][k]*dt + jrk[i][k]*dt*dt/2;
		}
	}
}

/*-----------------------------------------------------------------------------
 *  correct_step  --  takes one iteration to improve the new values of position
 *                    and velocities, effectively by using a higher-order
 *                    Taylor series constructed from the terms up to jrk at
 *                    the beginning and the end of the time step.
 *-----------------------------------------------------------------------------
 */

void correct_step(real pos[][NDIM], real vel[][NDIM],
				  const real acc[][NDIM], const real jrk[][NDIM],
				  const real old_pos[][NDIM], const real old_vel[][NDIM],
				  const real old_acc[][NDIM], const real old_jrk[][NDIM],
				  int n, real dt){
	for(int i = 0; i < n; i++){
		for(int k = 0; k < NDIM; k++){
			vel[i][k] = old_vel[i][k] + (old_acc[i][k] + acc[i][k])*dt/2
									  + (old_jrk[i][k] - jrk[i][k])*dt*dt/12;
			pos[i][k] = old_pos[i][k] + (old_vel[i][k] + vel[i][k])*dt/2
									  + (old_acc[i][k] - acc[i][k])*dt*dt/12;
		}
	}
}

/*-----------------------------------------------------------------------------
 *  get_acc_jrk_pot_coll  --  calculates accelerations and jerks, and as side
 *                            effects also calculates potential energy and
 *                            the time scale for significant changes in
 *                            local configurations to occur.
 *                                              __                          __
 *               GM_j         ;                |           r_ji . v_ji        |
 *    a_ji ==  --------_  rji ;  j_ji ==  GM_j | v_ji - 3 -------------  r_ji |
 *             |r_ji|^3       ;                |__           |r_ji|^3       __|
 *
 *  note: it would be cleaner to calculate potential energy and collision time
 *        in a separate function.  However, the current function is by far the
 *        most time consuming part of the whole program, with a double loop
 *        over all particles that is executed every time step.  Splitting off
 *        some of the work to another function would significantly increase
 *        the total computer time (by an amount close to a factor two).
 *
 *  We visit all particle pairs in a double {i,j} loop.  Acceleration,
 *  jerk, and potential energy, are calculated by adding successive terms;
 *  the estimate for the collision time is found by determining the minimum
 *  value over all particle pairs and over the two choices of collision time,
 *  position/velocity and sqrt(position/acceleration).
 *-----------------------------------------------------------------------------
 */

void get_acc_jrk_pot_coll(const real mass[], const real pos[][NDIM],
						  const real vel[][NDIM], real acc[][NDIM],
						  real jrk[][NDIM], real dst[], int n,
						  real & epot, real & coll_time){

	epot = 0;                         // potential energy
	real coll_time_q = DBL_MAX;       // collision time estimate to 4th power (quartic)

	for(int i = 0; i < n; i++){
		for(int k = 0; k < NDIM; k++){
			acc[i][k] = jrk[i][k] = 0;
		}
	}

	int p = 0;
	for(int i = 0; i < n; i++){
		for(int j = i+1; j < n; j++, p++){
			real rji[NDIM];           // vector from particle i to particle j
			real vji[NDIM];           // vji = d rji / d t

			real r2 = 0;              // | rji |^2
			real v2 = 0;              // | vji |^2
			real rv_r2 = 0;           // ( rij . vij ) / | rji |^2

			for(int k = 0; k < NDIM; k++){
				rji[k] = pos[j][k] - pos[i][k];
				vji[k] = vel[j][k] - vel[i][k];

				r2 += rji[k] * rji[k];
				v2 += vji[k] * vji[k];
				rv_r2 += rji[k] * vji[k];
			}

			rv_r2 /= r2;
			real r = sqrt(r2);        // | rji |
			real r3 = r * r2;         // | rji |^3

			dst[p] = r;  // record pairwise distances for output

			// add the {j (i)} contribution to the {i (j)} values of acceleration and jrk:

			real da2 = 0;  // square of the pairwise acceleration, used in collision estimates
			for(int k = 0; k < NDIM; k++){
				real da = rji[k] / r3;                           // see equations
				real dj = (vji[k] - 3 * rv_r2 * rji[k]) / r3;    // in the header

				da2 += da*da;

				acc[i][k] += mass[j] * da;  // using symmetry,
				acc[j][k] -= mass[i] * da;  // find pairwise
				jrk[i][k] += mass[j] * dj;  // acceleration
				jrk[j][k] -= mass[i] * dj;  // and jerk
			}

			// add the {i,j} contribution to the total potential energy for the system:
			epot -= mass[i] * mass[j] / r;

			// Estimate collision time by two methods, and take the smallest

			real coll_est_q;

			// first collision time estimate, based on unaccelerated linear motion:
			coll_est_q = (r2*r2) / (v2*v2);
			if(coll_time_q > coll_est_q){
				coll_time_q = coll_est_q;
			}

			// second collision time estimate, based on free fall:
			double mij = mass[i] + mass[j];  // add mass factors to finish calculating pairwise acceleration
			coll_est_q = G*r2/(da2*mij*mij); // masses already include a factor of G, so we cancel out the extra
			if(coll_time_q > coll_est_q){
				coll_time_q = coll_est_q;
			}
		}
	}                                     // from q for quartic back
	coll_time = sqrt(sqrt(coll_time_q));  // to linear collision time
}