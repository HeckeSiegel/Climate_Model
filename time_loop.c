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

int main()
{
    double t = 0.;
    double p0 = 1000.0; //hPa
    int nlev = 11;
    int nlyr = nlev - 1;
    int nwvl = 3;
    double tau_total = 1.;
    int dt = 15*60;

    double p[nlev],plyr[nlyr],z[nlev];
    double Tnlev[nlev],T[nlyr],theta[nlyr],tmp_T[nlyr];
    double Eup[nlev],Edown[nlev],deltaE[nlyr],tmp_deltaE[nlyr];

    //pressure profile
    for (int i = 0; i < nlev; i++) {
        p[i] = p0 * (double) i / (double) nlyr;
    }
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
    while(1==1){
        // gray atmosphere
        //oneBandAtmosphere(T,nlev,nlyr,tau_total,Edown,Eup);
        //window atmosphere
	threeBandAtmosphere(nwvl,nlyr,nlev,T,tau_total,Edown,Eup);
	//delta E
    	dE(deltaE,Edown,Eup,nlyr);
	for(int inlyr = 0; inlyr<nlyr; inlyr++){
	    tmp_T[inlyr] = T[inlyr];
	}
    	//deltaT
    	for (int i = 0; i < nlyr; i++){
	    T[i] += E2T(deltaE[i],100.,dt);
	    theta[i] = T2theta(T[i], plyr[i]);
    	}
	//convection
        convection(theta,nlyr);
        for (int i=0; i < nlyr; i++){
            T[i] = theta2T(theta[i],plyr[i]);
        }
	//time step
          //calculate absolute values
	for(int inlyr=0; inlyr<nlyr; inlyr++){
	    tmp_deltaE[inlyr] = tmp_T[inlyr]-T[inlyr];
	    if(tmp_deltaE[inlyr]<0.){
	        tmp_deltaE[inlyr] = -tmp_deltaE[inlyr];  
	    }
	}
          //find biggest change
	convection(tmp_deltaE,nlyr);
	while(tmp_deltaE[0] < 0.5){
	    if(dt > 60*60*24){
	        break;
	    }
	    dt += 15*60;
	}
	if(tmp_deltaE[0]<0.0001){
	    break;
	}
	else{
	    t+=dt;
	}
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
        printf("%6.3f %6.6f %6.6f\n",z[i]*1e-3,T[i],theta[i]);
    }
}


