#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "thermal_radiation.h"
#include <unistd.h>
#include "gnuplot_i.h"
#include "ascii.h"

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
//################################## main time loop that is being used for everything else #################################################
double timeLoop(double *T, double *theta, double *z, int nlyr, int nlev, int nwvl, double tau_total, double *Edown, double *Eup, double *wvl, double **tau){
    double t = 0.;
    double p0 = 1000.0; //hPa
    double equi = 1e-4; //threshold for breaking the time loop

    double p[nlev],plyr[nlyr];
    double Tnlev[nlev],deltaT[nlyr],dt,Ttoa;
    double deltaE[nlyr];

    FILE *file;
    int k = 0;

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
	Ttoa = T[nlyr-1]; // to compare temperature between 2 time steps
	char buffer[64]; // filename buffer
	snprintf(buffer, sizeof(char) * 64, "./plot/3band_animation/ascii/grey_animation_%i.asc", k);
	
        // gray atmosphere
        //oneBandAtmosphere(T,nlev,nlyr,tau_total,Edown,Eup,deltaE);
        //window atmosphere
	threeBandAtmosphere(nwvl,nlyr,nlev,T,tau_total,Edown,Eup,deltaE);
        // line by line atmosphere
        //multiBandAtmosphere(wvl,nwvl,nlyr,nlev,T,tau,Edown,Eup,deltaE);

	//################ find dt ####################################
	for(int inlyr=0; inlyr<nlyr; inlyr++){
		deltaT[inlyr] = E2T(deltaE[inlyr],dp);
		if(deltaT[inlyr]<0.){deltaT[inlyr] = -deltaT[inlyr];}
	}
	convection(deltaT,nlyr);
	dt = 1./deltaT[0];
	if(dt > 60*60*12){dt = 60*60*12;}
	//################ new T with calculated dt ####################
	for(int inlyr=0; inlyr<nlyr; inlyr++){
		T[inlyr] += E2T(deltaE[inlyr],dp)*dt;
	}
	//################# convection #################################
    	for (int inlyr = 0; inlyr < nlyr; inlyr++){
	    theta[inlyr] = T2theta(T[inlyr], plyr[inlyr]);
    	}
    	convection(theta,nlyr);
    	for (int inlyr=0; inlyr < nlyr; inlyr++){
    	    T[inlyr] = theta2T(theta[inlyr],plyr[inlyr]);
    	}
	t += dt;
	//p to z
        z[nlev-1] = 0.;
        z[nlev-2] = dp2dz(p[nlev-1]-plyr[nlyr-1],plyr[nlyr-1],T[nlyr-1]);
        for(int i=nlev-3; i>=0; i--){
	    z[i] = z[i+1] + dp2dz(plyr[i+1]-plyr[i],plyr[i],T[i]);
        }
	//print into file
	file = fopen(buffer, "wb");
	for (int i=0; i < nlyr; i++){
            fprintf(file,"%6.3f %6.3f\n",z[i]*1e-3,T[i]);
        }
	fclose(file );
	k++;
	//if the change in toa temperature is small enough the system is in equilibrium
	if(((Ttoa-T[nlyr-1]) < equi && (Ttoa-T[nlyr-1]) > 0) || ((Ttoa-T[nlyr-1]) > -equi && (Ttoa-T[nlyr-1]) < 0)){
	    break;
	}
	
    }
    return t;
}
//######################## plots height against temperature for every layer in every time step #######################
double plotTimeLoop(double *T, double *theta, double *z, int nlyr, int nlev, int nwvl, double tau_total, double *Edown, double *Eup, double *wvl, double **tau){
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
	Tsurf = T[0];
        // gray atmosphere
        //oneBandAtmosphere(T,nlev,nlyr,tau_total,Edown,Eup,deltaE);
        //window atmosphere
	//threeBandAtmosphere(nwvl,nlyr,nlev,T,tau_total,Edown,Eup,deltaE);
        // line by line atmosphere
        multiBandAtmosphere(wvl,nwvl,nlyr,nlev,T,tau,Edown,Eup,deltaE);

	//find dt
	for(int inlyr=0; inlyr<nlyr; inlyr++){
		deltaT[inlyr] = E2T(deltaE[inlyr],dp);
		if(deltaT[inlyr]<0.){deltaT[inlyr] = -deltaT[inlyr];}
	}
	convection(deltaT,nlyr);
	dt = 1./deltaT[0];
	if(dt > 60*60*12){dt = 60*60*12;}
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
      		gnuplot_setstyle   (g1, "linespoints");     /* draw lines and points */
      		gnuplot_set_xlabel (g1, "temperature [K]"); /* xaxis label */
      		gnuplot_set_ylabel (g1, "altitude [m]");    /* yaxis label */
      
      		/* plot temperature T as function of z and label with temperature */
      		gnuplot_plot_xy   (g1, T, z, nlyr, "Temperature") ;
		sleep(1);
    	 }
	t += dt;
	printf("T(surface) = %6.6f , z(surface in km) = %6.6f\n", T[nlyr-1],z[nlyr-1]*1e-3);
	//if the change in surface temperature is small enough the system is in equilibrium
	if(((Tsurf-T[0]) < equi && (Tsurf-T[0]) > 0) || ((Tsurf-T[0]) > -equi && (Tsurf-T[0]) < 0)){
	    break;
	}
    }
    /* close plot */
    sleep(50);
    gnuplot_close (g1) ;
    return t;
}
//#################################################################calculating the time loop for different taus###############################################
double plotT2tau(double *T, double *theta, double *z, int nlyr, int nlev, int nwvl, double tau_min, double tau_max, int tau_size, double *Edown, double *Eup, double *wvl, double **taulbl){
    double t = 0.;
    double tau_step = (tau_max - tau_min )/(tau_size-1);
    double tau[tau_size],Tequ[tau_size];
    
    gnuplot_ctrl *g1;
    g1 = gnuplot_init();

    int count = 0;
    for (double tau_total=tau_min; tau_total<=tau_max; tau_total += tau_step){
        t = timeLoop(T,theta,z,nlyr,nlev,nwvl,tau_total,Edown,Eup,wvl,taulbl);
        tau[count] = tau_total;
        Tequ[count] = T[nlyr-1];
        printf("Tsurf(%6.2f) = %6.3f\n",tau[count],Tequ[count]);
	count++;
    }
    //gnuplot_resetplot  (g1);  /* start with new plot rather than plotting into exisiting one */
    gnuplot_setstyle   (g1, "linespoints");      /* draw lines and points */
    gnuplot_set_xlabel (g1, "optical thickness tau");  /* xaxis label */
    gnuplot_set_ylabel (g1, "surface temperature [K]");    /* yaxis label */
    /* plot temperature T as function of z and label with temperature */
    gnuplot_plot_xy   (g1, tau, Tequ, count, "Temperature") ;
    
    sleep(10);	 
    /* close plot */
    //gnuplot_close (g1) ;

    return t;
}
int main()
{ 
    //for grey and 3-band atmospheres
    int nlyr = 50;
    int nlev = nlyr + 1;
    int nwvl = 3;
    double tau_total = 1.8;
    double z[nlev];
    double T[nlyr],theta[nlyr],t;
    double Eup[nlev],Edown[nlev];

    //for plotting different taus in grey and 3-band atmospheres
    double tau_min = 0.1;
    double tau_max = 20.1;
    int tau_size = (tau_max-tau_min)/1. + 1;
    
    //for line by line atmosphere
    int nwvlco2 = 0;
    int nwvlh2o = 0;
    int nyco2 = 0;
    int nyh2o = 0;
    double *wvnco2 = NULL; 
    double **tauco2 = NULL; 
    double *wvnh2o = NULL; 
    double **tauh2o = NULL;

    ASCII_file2xy2D ("lbl.arts/lbl.co2.asc", &nwvlco2, &nyco2, &wvnco2, &tauco2);
    ASCII_file2xy2D ("lbl.arts/lbl.h2o.asc", &nwvlh2o, &nyh2o, &wvnh2o, &tauh2o);

    int nlyrlbl = nyco2;
    int nlevlbl = nlyrlbl + 1;
    double zlbl[nlevlbl];
    double Tlbl[nlyrlbl],thetalbl[nlyrlbl];
    double Euplbl[nlevlbl],Edownlbl[nlevlbl];
    //add up tau's from h2o and co2
    for (int inwvl = 0; inwvl<nwvlco2; inwvl++){
	for (int inlyr = 0; inlyr<nlyrlbl; inlyr++){
		tauco2[inwvl][inlyr] = 0.*tauco2[inwvl][inlyr] + tauh2o[inwvl][inlyr];
	}
    }
    //############################ call different functions ##########################################

    t = timeLoop(T,theta,z,nlyr,nlev,nwvl,tau_total,Edown,Eup,wvnco2,tauco2); //for grey and window atmosphere
    //t = timeLoop(Tlbl,thetalbl,zlbl,nlyrlbl,nlevlbl,nwvlco2,tau_total,Edownlbl,Euplbl,wvnco2,tauco2); //for line by line atmosphere
    //t = plotTimeLoop(T,theta,z,nlyr,nlev,nwvl,tau_total,Edown,Eup,wvnco2,tauco2);
    //t = plotTimeLoop(Tlbl,thetalbl,zlbl,nlyrlbl,nlevlbl,nwvlco2,tau_total,Edownlbl,Euplbl,wvnco2,tauco2);
    //t = plotT2tau(T,theta,z,nlyr,nlev,nwvl,tau_min,tau_max,tau_size,Edown,Eup);

    //################################################################################################
    
    // print the output into a file
    /*FILE *fp;
    fp = fopen("test_lbl.txt", "w+");
    fprintf(fp, "after %6.2f years\n",t/(365.*24*60*60)); //time until equilibrium
    fprintf(fp,"Eup[TOA] = %6.3f\n",Eup[0]);             //Eup at top of atmosphere, should be Eearth when equilibrium is reached
    fprintf(fp,"z[km], T[K]:\n");                        //height and temperature profile
    for (int i=0; i < nlyr; i++){
        fprintf(fp,"%6.3f %6.3f %6.3f\n",z[i]*1e-3,T[i],theta[i]);
    }
    fclose(fp);*/
    ASCII_free_double(tauh2o, nwvlco2);
    ASCII_free_double(tauco2, nwvlh2o);
    free(wvnco2);
    free(wvnh2o);
}


