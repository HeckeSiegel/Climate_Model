#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "thermal_radiation.h"

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
    int nlev = 11;
    int nlyr = nlev - 1;
    int nwvl = 3;

    double B_layer[nlyr],B_surface,T[nlyr],tau[nlyr],tau0[nlyr],Eup[nlev],Edown[nlev],tmp_Eup[nlev],tmp_Edown[nlev];
    
    //second tau band
    for(int inlyr=0; inlyr<nlyr; inlyr++){
	tau0[inlyr] = 0.;
    }
    //wavelength band
    double wvlband[3][2] = {{1e-6,8e-6},{8e-6,12e-6},{12e-6,1e-2}};
    //test temperature profile
    double Tnlev[11] = {127.28,187.09,212.42,229.22,242.03,252.48,261.37,269.13,276.04,282.29,288.00};
    for (int i=0; i < nlyr; i++) {
        T[i] = (Tnlev[i]+Tnlev[i+1])/2.0;
    }
    printf("irradiance as a function of optical thickness of 1st and 3rd band\n");
    printf("tau, E_dn(SUR)[W/m2],E_up(TOA)[W/m2]:\n");
    for(double dtau=0.; dtau<4.8; dtau+=0.1){
        for (int inlyr = 0; inlyr< nlyr; inlyr++){
            tau[inlyr] = dtau/10.;
        }
        for (int inlev=0; inlev<nlev; inlev++){
            Eup[inlev] = 0.;
            Edown[inlev] = 0.;
        }
        for (int iwvl=0; iwvl<nwvl; iwvl++){
            B_surface = B_int(wvlband[iwvl][0],wvlband[iwvl][1],Tnlev[nlev-1]);
            for (int inlyr = 0; inlyr< nlyr; inlyr++){
                B_layer[inlyr] = B_int(wvlband[iwvl][0],wvlband[iwvl][1],T[inlyr]);
            }
	    if(iwvl==1){ 
	        schwarzschild(nlev,tau0,B_layer,B_surface,tmp_Edown,tmp_Eup);
            }
	    else{
    	        schwarzschild(nlev,tau,B_layer,B_surface,tmp_Edown,tmp_Eup);
	    }
	    for(int inlev=0; inlev<nlev; inlev++){
	        Eup[inlev] += tmp_Eup[inlev];
                Edown[inlev] += tmp_Edown[inlev];
	    }
        }
    
    printf("%6.1f %6.3f %6.3f \n",dtau,Edown[nlev-1],Eup[0]);
    }
}
