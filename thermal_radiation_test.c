#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "thermal_radiation.h"
#include "planck.h"

#define RA 287.0
#define CP 1004.0
#define g 9.81
#define Eearth 235. // W/m^2
#define sigma 5.67e-8
#define pi 3.14

double T2theta(double T, double p){
    double kappa = RA/CP;
    return T * pow(1000.0/p, kappa);
}
double theta2T(double theta, double p){
    double kappa = RA/CP;
    return theta * pow(p/1000.0 , kappa);
}
double E2T(double E, double p, double t){
    return E*g*t/(100*p*CP);
}
double dp2dz(double dp, double p, double T){
    return dp*T*RA/(g*p);
}
void convection(double *theta, int size){
    while(1==1){
        int stable = 1;
        for (int i  = size - 2; i >= 0; i--) {
            if (theta[i] < theta[i+1]) {
                double tmp = theta[i];
                theta[i] = theta[i+1];
                theta[i+1] = tmp;
                stable = 0;
            }
        }
        if(stable){
            break;
        }
        }
}

int main()
{
    double p0 = 1000.0; //hPa
    int nlev = 11;
    int nlyr = nlev - 1;
    int nwvl = 1;

    double p[nlev],plyr[nlyr],z[nlev],B_layer[nwvl][nlyr],B_surface[nwvl],T[nlyr],tau[nwvl][nlyr],Eup[nlev],Edown[nlev];
 
    //wavelength band
    double wvlband[3][2] = {{0.,8e-6},{8e-6,12e-6},{12e-6,130e-6}};
    //pressure profile
    for (int i = 0; i < nlev; i++) {
        p[i] = p0 * (double) i / (double) nlyr;
    }
    //tau profile
    int window = 0;
    double tau_total = 1.;
    double dtau = tau_total/nlyr;
    for (int iwvl=0; iwvl<nwvl; iwvl++){
        for (int inlyr = 0; inlyr< nlyr; inlyr++){
            tau[iwvl][inlyr] = dtau;
        }
    }
    if(window){
        for (int inlyr=0; inlyr<nlyr; inlyr++){
	    tau[1][inlyr] = 0.;
	}
    }
    //test temperature profile
    double Tnlev[11] = {127.28,187.09,212.42,229.22,242.03,252.48,261.37,269.13,276.04,282.29,288.00};
    for (int i=0; i < nlyr; i++) {
        T[i] = (Tnlev[i]+Tnlev[i+1])/2.0;
        plyr[i]=(p[i] + p[i+1]) / 2.0;
    }
    //Planck layer profile for each wavelength band
    // nwvl=1 means gray atmosphere
    if(nwvl == 1){
        B_surface[0] = B_gray(Tnlev[nlev-1]);
        for(int inlyr=0; inlyr<nlyr; inlyr++){
	    B_layer[0][inlyr] = B_gray(T[inlyr]);
	}
    }
    else{
        for (int iwvl=0; iwvl<nwvl; iwvl++){
            B_surface[iwvl] = B_int(wvlband[iwvl][0],wvlband[iwvl][1],Tnlev[nlev-1]);
        }
        for (int inlyr = 0; inlyr< nlyr; inlyr++){
	    for (int iwvl=0; iwvl<nwvl; iwvl++){
                B_layer[iwvl][inlyr] = B_int(wvlband[iwvl][0],wvlband[iwvl][1],T[inlyr]);
	    }
        }
    }
    //p to z
    z[nlev-1] = 0.;
    for(int i=nlev-2; i>=0; i--){
	z[i] = z[i+1] + dp2dz(100,plyr[i],T[i]);
    }
     //z, p, T
    printf("z[km], P[hPa], T[K]:\n");
    for (int i=0; i < nlev; i++){
        printf("%6.3f %6.1f %6.2f\n",z[i]*1e-3,p[i],Tnlev[i]);
    }
    //test
    double tmpEup[nlev], tmpEdown[nlev], tmpB_surface, tmpB_layer[nlyr], tmp_tau[nlyr];
    for(int inlev=0; inlev<nlev; inlev++){
	Eup[inlev] = 0.;
	Edown[inlev] = 0.;
    }
    for(int iwvl=0; iwvl<nwvl; iwvl++){
	for(int inlyr=0; inlyr<nlyr; inlyr++){
	    tmpB_layer[inlyr] = B_layer[iwvl][inlyr];
	    tmp_tau[inlyr] = tau[iwvl][inlyr];
	}
	tmpB_surface = B_surface[iwvl];
    	schwarzschild(nlev,tmp_tau,tmpB_layer,tmpB_surface,tmpEdown,tmpEup);
	for(int inlev=0; inlev<nlev; inlev++){
	    Eup[inlev] += tmpEup[inlev];
	    Edown[inlev] += tmpEdown[inlev];
	}
    }
    printf("tau_total = %2f :\n", tau_total);
    printf("z[km], P[hPa], E_dn[W/m2], E_up[W/m2]:\n");
    for (int i=0; i < nlev; i++){
        printf("%6.3f %6.1f %6.3f %6.3f \n",z[i]*1e-3,p[i],Edown[i],Eup[i]);
    }
}
