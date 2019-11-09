#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "thermal_radiation.h"

#define RA 287.0
#define CP 1005.0
#define g 9.81
#define Eearth 240. // W/m^2
#define sigma 5.670374e-8

double T2theta(double T, double p){
    double kappa = RA/CP;
    return T * pow(1000.0/p, kappa);
}
double theta2T(double theta, double p){
    double kappa = RA/CP;
    return theta * pow(p/1000.0 , kappa);
}
double E2T(double E, double dp, double t){
    return E*g*t/(100*dp*CP);
}
double dp2dz(double dp, double p, double T){
    return dp*T*RA/(g*p);
}
double boltzmann(double T){
    return sigma*pow(T,4);
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
void dT(double *deltaE, double *T, int nlyr ,double dp, double dt){
    for (int inlyr = 0; inlyr < nlyr; inlyr++){
	T[inlyr] += E2T(deltaE[inlyr],dp,dt);
    }
}
void absDeltaT(double *T, double *tmp_T, double *deltaT, int nlyr){
    for(int inlyr=0; inlyr<nlyr; inlyr++){
	deltaT[inlyr] = T[inlyr] - tmp_T[inlyr];
	if(deltaT[inlyr]<0.){
		deltaT[inlyr] = -deltaT[inlyr];
	}
    }
}
int main()
{
    double t = 0.;
    int years = 2;
    double p0 = 1000.0; //hPa
    int nlev = 11;
    int nlyr = nlev - 1;
    int nwvl = 3;
    double tau_total = 1.0;


    double p[nlev],plyr[nlyr],z[nlev];
    double Tnlev[nlev],T[nlyr],theta[nlyr],tmp_T[nlyr],deltaT[nlyr],dt;
    double Eup[nlev],Edown[nlev],deltaE[nlyr];

    //pressure profile
    for (int i = 0; i < nlev; i++) {
        p[i] = p0 * (double) i / (double) nlyr;
    }
    double dp = p[1]-p[0];
    //random level temperature profile
    for (int i=0; i < nlev; i++) {
        Tnlev[i] = 200.0 + ((double)random() / (double) RAND_MAX)*100;
    }
    //layer temperature, potential temperature, pressure
    for (int i=0; i < nlyr; i++) {
        T[i] = (Tnlev[i] + Tnlev[i+1]) / 2.0;
        plyr[i]=(p[i] + p[i+1]) / 2.0;
        theta[i] = T2theta(T[i], plyr[i]);
    }
    //time loop
    while(t<years*60*60*24*365){
	dt = 15*60;
	for(int inlyr=0; inlyr<nlyr; inlyr++){
	    tmp_T[inlyr] = T[inlyr];
	}
        // gray atmosphere
        //oneBandAtmosphere(T,nlev,nlyr,tau_total,Edown,Eup,deltaE);
        //window atmosphere
	threeBandAtmosphere(nwvl,nlyr,nlev,T,tau_total,Edown,Eup,deltaE);

	//find right time step
	while(1==1){
	    //new T
	    dT(deltaE,T,nlyr,dp,dt);
	    //convection
    	    for (int inlyr = 0; inlyr < nlyr; inlyr++){
		theta[inlyr] = T2theta(T[inlyr], plyr[inlyr]);
    	    }
    	    convection(T,nlyr);
    	    for (int inlyr=0; inlyr < nlyr; inlyr++){
        	T[inlyr] = theta2T(theta[inlyr],plyr[inlyr]);
    	    }
	    //absolute values of temperature difference
	    absDeltaT(T, tmp_T, deltaT, nlyr);
	    //find biggest temperature change
	    convection(deltaT,nlyr);
	    if(deltaT[0]>1.0){
		break;
	    }
	    else{
		dt += 15*60;
	    }
	    if(dt>60*60*24){break;}
	}
	t += dt;
	printf("%6.3f %6.2f %6.2f\n",T[nlyr-1],deltaE[nlyr-1],deltaE[0]);
    }
    //p to z
    z[nlev-1] = 0.;
    for(int i=nlev-2; i>=0; i--){
	z[i] = z[i+1] + dp2dz(100,plyr[i],T[i]);
    }
    printf("tau_total =%6.1f after %6.1f years\n",tau_total,t/(365.*24*60*60));
    printf("Eup[TOA] = %6.3f\n",Eup[0]);
    printf("z[km], T, theta:\n");
    for (int i=0; i < nlyr; i++){
        printf("%6.3f %6.3f %6.3f\n",z[i]*1e-3,T[i],theta[i]);
    }
}


