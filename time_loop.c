#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "thermal_radiation.h"

#define RA 287.0
#define CP 1004.0
#define g 9.81
#define Eearth 240.4 // W/m^2
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
    int t = 0;
    double p0 = 1000.0; //hPa
    int nlev = 11;
    int nlyr = nlev - 1;
    int years = 2;
    int nwvl = 3;
    double tau_total = 1.0;
    int dt = 15*60;

    double p[nlev],plyr[nlyr],z[nlev];
    double Tnlev[nlev],T[nlyr],theta[nlyr];
    double Eup[nlev],Edown[nlev],deltaE[nlyr];

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
    while(t<years*365*24*3600){
        /* gray atmosphere
        oneBandAtmosphere(T,nlev,nlyr,tau_total,Edown,Eup);
        */

        //window atmosphere
	threeBandAtmosphere(nwvl,nlyr,nlev,T,tau_total,Edown,Eup);
	
    	dE(deltaE,Edown,Eup,nlyr);

        //time step
	double gradE = deltaE[nlyr-4] - deltaE[nlyr-5];
        if(gradE<0){
	    gradE = -gradE;  
	}
	while(gradE*g*dt/(10000.*CP) < 1.0){
	    dt += 10*60;
	    if(dt > 60*60*24){
	        break;
	    }
	}
    	//deltaT
    	for (int i = 0; i < nlyr; i++){
	    T[i] += E2T(deltaE[i],p[i+1]-p[i],dt);
	    theta[i] = T2theta(T[i], plyr[i]);
    	}
	//heating
        T[nlyr-1] += E2T(Eearth,p[nlev-1]-p[nlev-2],dt);
	//convection
        convection(theta,nlyr);
        for (int i=0; i < nlyr; i++){
        T[i] = theta2T(theta[i],plyr[i]);
        }
	//printf("%2d\n",dt/60);
	t+=dt;
    }
    //p to z
    z[nlev-1] = 0.;
    for(int i=nlev-2; i>=0; i--){
	z[i] = z[i+1] + dp2dz(100,plyr[i],T[i]);
    }
    printf("tau_total =%6.1f after %2d years\n",tau_total,years);
    printf("Eup[TOA] = %6.3f\n",Eup[0]);
    printf("Edown[TOA] = %6.3f\n",Edown[0]);
    printf("z[km], T, theta:\n");
    for (int i=0; i < nlyr; i++){
        printf("%6.3f %6.6f %6.6f\n",z[i]*1e-3,T[i],theta[i]);
    }
}

