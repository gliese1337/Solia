#ifndef EVOLVE_H
#define EVOLVE_H

void evolve_step(const double mass[], double pos[][NDIM], double vel[][NDIM],
				 double acc[][NDIM], double jrk[][NDIM], double dst[],
				 double old_pos[][NDIM], double old_vel[][NDIM],
				 double old_acc[][NDIM], double old_jrk[][NDIM],
				 int n, double dt, double & epot, double & coll_time);

void predict_step(double pos[][NDIM], double vel[][NDIM],
				  const double acc[][NDIM], const double jrk[][NDIM],
				  int n, double dt);

void correct_step(double pos[][NDIM], double vel[][NDIM],
				  const double acc[][NDIM], const double jrk[][NDIM],
				  const double old_pos[][NDIM], const double old_vel[][NDIM],
				  const double old_acc[][NDIM], const double old_jrk[][NDIM],
				  int n, double dt);

void get_acc_jrk_pot_coll(const double mass[], const double pos[][NDIM],
						  const double vel[][NDIM], double acc[][NDIM],
						  double jrk[][NDIM], double dst[], int n, double & epot,
						  double & coll_time);

#endif