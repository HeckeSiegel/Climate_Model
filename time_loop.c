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
#define h 6.626e-34 //Js
#define c 2.99e8 // m/s
#define kB 1.38e-23 // m^2kgs^-2K^-1

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
    int t = 0;
    int dt = 15*60; //in s, nicht 60 statt 16?

    double p0 = 1000.0; //hPa
    int nlev = 11;
    int nlyr = nlev - 1;
    int years = 1;
    int nwvl = 3;

    double p[nlev],z[nlev],plyr[nlyr],B_layer[nwvl][nlyr],B_surface[nwvl],Tnlev[nlev],T[nlyr],theta[nlyr],tau[nwvl][nlyr],Eup[nlyr],Edown[nlev],deltaE[nlyr];
    
    //wavelength band
    double wvlband[3][2] = {{0.,8e-6},{8e-6,12e-6},{12e-6,120e-6}};
    //pressure profile
    for (int i = 0; i < nlev; i++) {
        p[i] = p0 * (double) i / (double) nlyr;
    }
    //tau profile
    int window = 1;
    double tau_total = 1;
    double dtau = tau_total/nlyr;
    for (int iwvl=0; iwvl<wvl; iwvl++){
        for (int inlyr = 0; inlyr< nlyr; inlyr++){
            tau[iwvl][inlyr] = dtau;
        }
    }
    if(window){
        for (int inlyr=0; inlyr<nlyr; inylr++){
	    tau[1][inlyr] = 0;
	}
    {
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
    //p to z
    z[nlev-1] = 0.;
    for(int i=nlev-2; i>=0; i--){
	z[i] = z[i+1] + dp2dz(100,plyr[i],T[i]);
    }
    //time loop
    while(t<years*365*24*3600){
        t+=dt;
        //Planck layer profile for each wavelength band
	for (int iwvl=0; iwvl<wvl; iwvl++){
            B_surface[iwvl] = B_int(wvlband[iwvl][0],wvlband[iwvl][1],T[nlyr-1]);
	}
        for (int inlyr = 0; inlyr< nlyr; inlyr++){
	    for (int iwvl=0; iwvl<wvl; iwvl++){
                B_layer[iwvl][inlyr] = B_int(wvlband[iwvl][0],wvlband[iwvl][1],T[inlyr]);
	    }
        }
	//heating
        T[nlyr-1] += E2T(Eearth,plyr[nlyr-1],dt); //heating
        for(int i = 0; i<nlyr ; i++){
            theta[i] = T2theta(T[i], plyr[i]);
	}
	//radiation transport for each wavelength band TODO:run schwarzschild for each band, add up
    	    schwarzschild(nlev,tau,B_layer,B_surface,Edown,Eup);

    	dE(deltaE,Edown,Eup,nlyr);
    	//deltaT
    	for (int i = 0; i < nlyr; i++){
	    T[i] += E2T(deltaE[i],plyr[i],dt);
	    theta[i] = T2theta(T[i], plyr[i]);
    	}
	//convection
        convection(theta,nlyr);
        for (int i=0; i < nlyr; i++){
        T[i] = theta2T(theta[i],plyr[i]);
        }
    }
    printf("z[km], T, theta after %2d years:\n", years);
    for (int i=0; i < nlyr; i++){
        printf("%6.3f %6.1f %6.1f\n",z[i]*1e-3,T[i],theta[i]);
    }
}
