#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "thermal_radiation.h"
#include <unistd.h>
#include "gnuplot_i.h"

#define RA 287.0
#define CP 1005.0
#define g 9.81
#define Eearth 235. // W/m^2
#define sigma 5.670374e-8

double T2theta(double T, double p){
    double kappa = RA/CP;
    return T * pow(1000.0/p, kappa);
}
double theta2T(double theta, double p){
    double kappa = RA/CP;
    return theta * pow(p/1000.0 , kappa);
}
double E2T(double E, double dp){
    return E*g/(100*dp*CP);
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
double plotZ2T(double *T, double *theta, double *z, int nlyr, int nlev, int nwvl, double tau_total, double *Edown, double *Eup){
    double t = 0.;
    double p0 = 1000.0; //hPa
    double equi = 1e-5;

    double p[nlev],plyr[nlyr];
    double Tnlev[nlev],deltaT[nlyr],dt,Tsurf;
    double deltaE[nlyr];
    
    int tcounter=0;
    gnuplot_ctrl *g1;
    g1 = gnuplot_init();

    //pressure profile
    for (int inlev = 0; inlev < nlev; inlev++) {
        p[inlev] = p0 * (double) inlev / (double) nlyr;
    }
    double dp = p[1]-p[0];
    //random level temperature profile
    for (int inlev=0; inlev < nlev; inlev++) {
        Tnlev[inlev] = 200.0 + ((double)random() / (double) RAND_MAX)*100;
    }
    //layer temperature, pressure
    for (int inlev=0; inlev < nlyr; inlev++) {
        T[inlev] = (Tnlev[inlev] + Tnlev[inlev+1]) / 2.0;
        plyr[inlev]=(p[inlev] + p[inlev+1]) / 2.0;
    }
    //time loop
    while(1==1){
	Tsurf = T[nlyr-1];
	dt = 15*60;
        // gray atmosphere
        //oneBandAtmosphere(T,nlev,nlyr,tau_total,Edown,Eup,deltaE);
        //window atmosphere
	threeBandAtmosphere(nwvl,nlyr,nlev,T,tau_total,Edown,Eup,deltaE);
	//find dt
	for(int inlyr=0; inlyr<nlyr; inlyr++){
		deltaT[inlyr] = E2T(deltaE[inlyr],dp);
		if(deltaT[inlyr]<0.){deltaT[inlyr] = -deltaT[inlyr];}
	}
	convection(deltaT,nlyr);
	while(deltaT[0]*dt <= 1.){
		dt += 15*60;
		if(dt > 60*60*24){break;};
	}

	//new T with calculated dt
	for(int inlyr=0; inlyr<nlyr; inlyr++){
		T[inlyr] += E2T(deltaE[inlyr],dp)*dt;
	}
	//convection
    	for (int inlyr = 0; inlyr < nlyr; inlyr++){
	    theta[inlyr] = T2theta(T[inlyr], plyr[inlyr]);
    	}
    	convection(theta,nlyr);
    	for (int inlyr=0; inlyr < nlyr; inlyr++){
    	    T[inlyr] = theta2T(theta[inlyr],plyr[inlyr]);
    	}
	//p to z
    	z[nlev-1] = 0.;
	z[nlev-2] = dp2dz(p[nlev-1]-plyr[nlyr-1],plyr[nlyr-1],T[nlyr-1]);
    	for(int i=nlev-3; i>=0; i--){
		z[i] = z[i+1] + dp2dz(plyr[i+1]-plyr[i],plyr[i],T[i]);
    	}
	
	/* number of time steps */
        tcounter++;
	if (tcounter%10 == 0) {
      
      		gnuplot_resetplot  (g1);  /* start with new plot rather than plotting into exisiting one */
      		gnuplot_setstyle   (g1, "linespoints");      /* draw lines and points */
      		gnuplot_set_xlabel (g1, "temperature [K]");  /* xaxis label */
      		gnuplot_set_ylabel (g1, "altitude [m]");    /* yaxis label */
      
      		/* plot temperature T as function of z and label with temperature */
      		gnuplot_plot_xy   (g1, T, z, nlyr, "Temperature") ;
		//sleep(1);
    	 }
	t += dt;
	printf("T(surface) = %6.6f , z(surface in km) = %6.6f\n", T[nlyr-1],z[nlyr-1]*1e-3);
	//if the change in surface temperature is small enough the system is in equilibrium
	if(((Tsurf-T[nlyr-1]) < equi && (Tsurf-T[nlyr-1]) > 0) || ((Tsurf-T[nlyr-1]) > -equi && (Tsurf-T[nlyr-1]) < 0)){
	    break;
	}
    }
    /* close plot */
    gnuplot_close (g1) ;
    return t;
}
double timeLoop(double *T, double *theta, double *z, int nlyr, int nlev, int nwvl, double tau_total, double *Edown, double *Eup){
    double t = 0.;
    double p0 = 1000.0; //hPa
    double equi = 1e-5;

    double p[nlev],plyr[nlyr];
    double Tnlev[nlev],deltaT[nlyr],dt,Tsurf;
    double deltaE[nlyr];

    //pressure profile
    for (int inlev = 0; inlev < nlev; inlev++) {
        p[inlev] = p0 * (double) inlev / (double) nlyr;
    }
    double dp = p[1]-p[0];
    //random level temperature profile
    for (int inlev=0; inlev < nlev; inlev++) {
        Tnlev[inlev] = 200.0 + ((double)random() / (double) RAND_MAX)*100;
    }
    //layer temperature, pressure
    for (int inlev=0; inlev < nlyr; inlev++) {
        T[inlev] = (Tnlev[inlev] + Tnlev[inlev+1]) / 2.0;
        plyr[inlev]=(p[inlev] + p[inlev+1]) / 2.0;
    }
    //time loop
    while(1==1){
	Tsurf = T[nlyr-1];
	dt = 15.*60.;
        // gray atmosphere
        //oneBandAtmosphere(T,nlev,nlyr,tau_total,Edown,Eup,deltaE);
        //window atmosphere
	threeBandAtmosphere(nwvl,nlyr,nlev,T,tau_total,Edown,Eup,deltaE);

	//find dt
	for(int inlyr=0; inlyr<nlyr; inlyr++){
		deltaT[inlyr] = E2T(deltaE[inlyr],dp);
		if(deltaT[inlyr]<0.){deltaT[inlyr] = -deltaT[inlyr];}
	}
	convection(deltaT,nlyr);
	while(deltaT[0]*dt <= 1.){
		dt += 15*60;
		if(dt > 60*60*24){break;};
	}

	//new T with calculated dt
	for(int inlyr=0; inlyr<nlyr; inlyr++){
		T[inlyr] += E2T(deltaE[inlyr],dp)*dt;
	}
	//convection
    	for (int inlyr = 0; inlyr < nlyr; inlyr++){
	    theta[inlyr] = T2theta(T[inlyr], plyr[inlyr]);
    	}
    	convection(theta,nlyr);
    	for (int inlyr=0; inlyr < nlyr; inlyr++){
    	    T[inlyr] = theta2T(theta[inlyr],plyr[inlyr]);
    	}
	t += dt;
	printf("T(surface) = %6.6f , t(days) = %6.6f\n", T[nlyr-1],t/(60.*60.*24.));
	//if the change in surface temperature is small enough the system is in equilibrium
	if(((Tsurf-T[nlyr-1]) < equi && (Tsurf-T[nlyr-1]) > 0) || ((Tsurf-T[nlyr-1]) > -equi && (Tsurf-T[nlyr-1]) < 0)){
	    break;
	}
	
    }
    //p to z
    z[nlev-1] = 0.;
    z[nlev-2] = dp2dz(p[nlev-1]-plyr[nlyr-1],plyr[nlyr-1],T[nlyr-1]);
    for(int i=nlev-3; i>=0; i--){
	z[i] = z[i+1] + dp2dz(plyr[i+1]-plyr[i],plyr[i],T[i]);
    }
    return t;
}
double plotT2tau(double *T, double *theta, double *z, int nlyr, int nlev, int nwvl, double tau_min, double tau_max, int tau_size, double *Edown, double *Eup){
    double p0 = 1000.0; //hPa
    double equi = 1e-4;
    double t = 0;
    double tau_step = (tau_max - tau_min )/(tau_size-1);
    double p[nlev],plyr[nlyr];
    double Tnlev[nlev],deltaT[nlyr],dt,Tsurf,tau[tau_size],Tequ[tau_size];
    double deltaE[nlyr];
    
    gnuplot_ctrl *g1;
    g1 = gnuplot_init();

    //pressure profile
    for (int inlev = 0; inlev < nlev; inlev++) {
        p[inlev] = p0 * (double) inlev / (double) nlyr;
    }
    double dp = p[1]-p[0];
    //random level temperature profile
    for (int inlev=0; inlev < nlev; inlev++) {
        Tnlev[inlev] = 200.0 + ((double)random() / (double) RAND_MAX)*100;
    }
    //layer temperature, pressure
    for (int inlev=0; inlev < nlyr; inlev++) {
        T[inlev] = (Tnlev[inlev] + Tnlev[inlev+1]) / 2.0;
        plyr[inlev]=(p[inlev] + p[inlev+1]) / 2.0;
    }
    int count = 0;
    //time loop
    for (double tau_total=tau_min; tau_total<=tau_max; tau_total += tau_step){
    t = 0.;
    while(1==1){
	Tsurf = T[nlyr-1];
	dt = 15*60;
        // gray atmosphere
        //oneBandAtmosphere(T,nlev,nlyr,tau_total,Edown,Eup,deltaE);
        //window atmosphere
	threeBandAtmosphere(nwvl,nlyr,nlev,T,tau_total,Edown,Eup,deltaE);
	//find dt
	for(int inlyr=0; inlyr<nlyr; inlyr++){
		deltaT[inlyr] = E2T(deltaE[inlyr],dp);
		if(deltaT[inlyr]<0.){deltaT[inlyr] = -deltaT[inlyr];}
	}
	convection(deltaT,nlyr);
	while(deltaT[0]*dt <= 1.){
		dt += 15*60;
		if(dt > 60*60*24){break;};
	}

	//new T with calculated dt
	for(int inlyr=0; inlyr<nlyr; inlyr++){
		T[inlyr] += E2T(deltaE[inlyr],dp)*dt;
	}
	//convection
    	for (int inlyr = 0; inlyr < nlyr; inlyr++){
	    theta[inlyr] = T2theta(T[inlyr], plyr[inlyr]);
    	}
    	convection(theta,nlyr);
    	for (int inlyr=0; inlyr < nlyr; inlyr++){
    	    T[inlyr] = theta2T(theta[inlyr],plyr[inlyr]);
    	}
	//p to z
    	z[nlev-1] = 0.;
	z[nlev-2] = dp2dz(p[nlev-1]-plyr[nlyr-1],plyr[nlyr-1],T[nlyr-1]);
    	for(int i=nlev-3; i>=0; i--){
		z[i] = z[i+1] + dp2dz(plyr[i+1]-plyr[i],plyr[i],T[i]);
    	}
	t += dt;
        //printf("T(surface) = %6.6f , t(days) = %6.6f\n", T[nlyr-1],t/(60.*60.*24.));
	//if the change in surface temperature is small enough the system is in equilibrium
	if(((Tsurf-T[nlyr-1]) < equi && (Tsurf-T[nlyr-1]) > 0) || ((Tsurf-T[nlyr-1]) > -equi && (Tsurf-T[nlyr-1]) < 0)){
	    break;
	}
    }
    tau[count] = tau_total;
    Tequ[count] = T[nlyr-1];
    
    printf("Tsurf(%6.2f) = %6.3f\n",tau[count],Tequ[count]);
    
    }
    //gnuplot_resetplot  (g1);  /* start with new plot rather than plotting into exisiting one */
    gnuplot_setstyle   (g1, "linespoints");      /* draw lines and points */
    gnuplot_set_xlabel (g1, "optical thickness tau");  /* xaxis label */
    gnuplot_set_ylabel (g1, "surface temperature [K]");    /* yaxis label */
    /* plot temperature T as function of z and label with temperature */
    gnuplot_plot_xy   (g1, tau, Tequ, nlyr, "Temperature") ;
    sleep(10);	 
    /* close plot */
    gnuplot_close (g1) ;
    return t;
}
int main()
{
    int nlyr = 20;
    int nlev = nlyr + 1;
    int nwvl = 3;
    double tau_total = 98.0;
    double tau_min = 0.5;
    double tau_max = 2.0;
    double tau_size = 16;
    double z[nlev];
    double T[nlyr],theta[nlyr],t;
    double Eup[nlev],Edown[nlev];

    t = timeLoop(T,theta,z,nlyr,nlev,nwvl,tau_total,Edown,Eup);

    printf("tau_total = %6.1f after %6.2f years\n",tau_total,t/(365.*24*60*60));
    printf("Eup[TOA] = %6.3f\n",Eup[0]);
    printf("z[km], T[K]:\n");
    for (int i=0; i < nlyr; i++){
        printf("%6.3f %6.3f %6.3f\n",z[i]*1e-3,T[i],theta[i]);
    }
}


