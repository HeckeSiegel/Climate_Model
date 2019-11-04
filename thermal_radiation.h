#ifndef _thermal_radiation_h
#define _thermal_radiation_h 1

double B_gray(double T);
double B_int(double lambda_1, double lambda_2, double T);
double alpha(double mu, double tau);
void schwarzschild(int nlev, double *tau, double *B_layer, double B_surface, double *Edown, double *Eup);
void dE(double *deltaE, double *Edown, double *Eup, int nlyr);

#endif
