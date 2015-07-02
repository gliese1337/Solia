#ifndef NBODYIO_H
#define NBODYIO_H

void get_snapshot(double mass[], double pos[][NDIM], double vel[][NDIM], int n);

void put_snapshot(const double mass[], const double pos[][NDIM],
				  const double vel[][NDIM], const double dst[],
				  int n, double t);

void write_diagnostics(const double mass[], const double pos[][NDIM],
					   const double vel[][NDIM], const double acc[][NDIM],
					   const double jrk[][NDIM], int n, double t, double epot,
					   int nsteps, bool x_flag);

#endif