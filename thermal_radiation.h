#ifndef _thermal_radiation_h
#define _thermal_radiation_h 1

double B_gray(double T);
double B_int(double lambda_1, double lambda_2, double T);
double alpha(double mu, double tau);
void schwarzschild(int nlev, double *tau, double *B_layer, double B_surface, double *Edown, double *Eup);
void dE(double *deltaE, double *Edown, double *Eup, int nlyr);
void oneBandAtmosphere(double *T, int nlev, int nlyr, double tau_total, double *Edown, double *Eup, double *deltaE);
void threeBandAtmosphere(int nwvl, int nlyr, int nlev, double *T, double tau_total, double *Edown, double *Eup, double *deltaE);
void multiBandAtmosphere(double *wvl,int nwvl, int nlyr, int nlev, double *T, double **tau, double *Edown, double *Eup, double *deltaE);
void kAtmosphere(double **wgt_lw, double *band_lbound, double *band_ubound, int nbands, int nlyr, int nlev, double *T, double **tau, double *Edown, double *Eup, double *deltaE);
void eddington_v2 (double dtau, double g, double omega0, double mu0,
		   double *t, double *r, double *rdir, double *sdir, double *tdir);
void solar_rt (double *deltaE, double Ag, int nlev, int nlyr, double *wgt_sw, double *band_lbound, double *band_ubound, int nbands, double g[nbands][nlyr],double tau[nbands][nlyr],double omega0[nbands][nlyr]);
#endif
