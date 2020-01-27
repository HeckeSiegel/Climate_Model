#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "thermal_radiation.h"
#include <unistd.h>
#include "gnuplot_i.h"
#include "ascii.h"
#include "fpda_rrtm_lw.h"
#include "fpda_rrtm_sw.h"

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
double magnus(double T){
    return 6.112*exp(17.62*(T-273.15)/(243.12+T-273.15));
}
double deltaw(double w0, double g0){
    return (1-g0*g0)*w0/(1-w0*g0*g0);
}
double deltag(double g0){
    return g0/(1+g0*g0);
}
double deltatau(double tau, double w0, double g0){
    return tau*(1-w0*g0*g0);
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
//################################## time loop for grey, 3-band and line-by-line atmosphere ##############################
double timeLoop(double *T, double *theta, double *z, int nlyr, int nlev, int nwvl, double tau_total, double *Edown, double *Eup, double *wvl, double **tau){
    double t = 0.;
    double p0 = 1000.0; //hPa
    double equi = 1e-6; //threshold for breaking the time loop

    double p[nlev],plyr[nlyr];
    double Tnlev[nlev],deltaT[nlyr],dt,Ttoa;
    double deltaE[nlyr];

    //FILE *file;
    //int k = 0;

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
    while((t/(24*60*60)) < 260){
	Ttoa = T[nlyr-1]; // to compare temperature between 2 time steps

	//###################################### choose which kind of atmosphere to simulate #############################
        // gray atmosphere
        //oneBandAtmosphere(T,nlev,nlyr,tau_total,Edown,Eup,deltaE);
        //window atmosphere
	//threeBandAtmosphere(nwvl,nlyr,nlev,T,tau_total,Edown,Eup,deltaE);
        // line by line atmosphere
        multiBandAtmosphere(wvl,nwvl,nlyr,nlev,T,tau,Edown,Eup,deltaE);
        //################################################################################################################

	//################ find dt #######################################################################################
	for(int inlyr=0; inlyr<nlyr; inlyr++){
		deltaT[inlyr] = E2T(deltaE[inlyr],dp);
		if(deltaT[inlyr]<0.){deltaT[inlyr] = -deltaT[inlyr];}
	}
	convection(deltaT,nlyr);
	dt = 1./deltaT[0];
	if(dt > 60*60*12){dt = 60*60*12;}

	//################ new T with calculated dt #####################################################################
	for(int inlyr=0; inlyr<nlyr; inlyr++){
		T[inlyr] += E2T(deltaE[inlyr],dp)*dt;
	}
	T[nlyr-1] += E2T(Eearth,dp)*dt;
	//################# convection ##################################################################################
    	for (int inlyr = 0; inlyr < nlyr; inlyr++){
	    theta[inlyr] = T2theta(T[inlyr], plyr[inlyr]);
    	}
    	convection(theta,nlyr);
    	for (int inlyr=0; inlyr < nlyr; inlyr++){
    	    T[inlyr] = theta2T(theta[inlyr],plyr[inlyr]);
    	}
	t += dt;
	//################ convert p to z ##############################################################################
        z[nlev-1] = 0.;
        z[nlev-2] = dp2dz(p[nlev-1]-plyr[nlyr-1],plyr[nlyr-1],T[nlyr-1]);
        for(int i=nlev-3; i>=0; i--){
	    z[i] = z[i+1] + dp2dz(plyr[i+1]-plyr[i],plyr[i],T[i]);
        }
	printf("T(surf) = %6.3f\n", T[nlyr-1]);
	
	//################################### print into file to make plot ###############################################
	/*char buffer[64]; // filename buffer
	snprintf(buffer, sizeof(char) * 64, "./plot/3band_animation/ascii/grey_animation_%i.asc", k);*/
	//file = fopen("./plot/animation/....asc", "wb");
	//file = fopen("./plot/lbl_zero_all_new.asc", "ab");
        //fprintf(file,"%6.6f %6.6f\n",T[nlyr-1],t/(24*60*60));
	//fclose(file );
	//k++;
        //################################################################################################################

	//### if the change in toa temperature is small enough the system is in equilibrium ##############################
	/*if(((Ttoa-T[nlyr-1]) < equi && (Ttoa-T[nlyr-1]) > 0) || ((Ttoa-T[nlyr-1]) > -equi && (Ttoa-T[nlyr-1]) < 0)){
	    break;
	}*/
	
    }
    return t;
}
//################################## time loop for k distribution atmosphere ##############################################
double kTimeLoop(int nlyr, int nlev){
    double t = 0.;
    double p0 = 1000.0; //hPa
    double equi = 1e-5; //threshold for breaking the time loop
	
    //FILE *file;

    double p[nlev], plyr[nlyr];
    double T[nlyr], theta[nlyr], z[nlyr], deltaT[nlyr], Ttoa, dt;
    double Edown[nlev], Eup[nlev], deltaE_clear[nlyr], deltaE_cloud[nlyr];
    double deltaE_sw_clear[nlyr],deltaE_sw_cloud[nlyr];
    double h2ovmr[nlyr]   , o3vmr[nlyr]    , co2vmr[nlyr]   , ch4vmr[nlyr]   , n2ovmr[nlyr] , o2vmr[nlyr];
    double cfc11vmr[nlyr] , cfc12vmr[nlyr] , cfc22vmr[nlyr] , ccl4vmr[nlyr];
    
    //for longwave radiation
    int nbands;
    double * band_lbound;   // [nbands]    
    double * band_ubound;   // [nbands]    
    double **tau;   // [nbands][nlyr]       
    double **wgt_lw;        // [nbands][nlyr]       
    
    // for shortwave radiation
    int nbands_sw;
    double *band_lbound_sw; // [nbands]    
    double *band_ubound_sw; // [nbands]    

    double **dtau_mol_sw;   // [nbands][nlay]       
    double **dtau_ray_sw;   // [nbands][nlay]       
    double *wgt_sw;         // [nbands]
    double g0 = 0.;
    double Ag = 0.3;
    
    //read in h2o and O3 concentrations
    int nrows = 0;
    double *tmp1 = NULL;
    double *tmp2 = NULL; 
    double *tmp3 = NULL;  
    double *h20 = NULL;
    double *o3 = NULL;
    double rel_F[nlev];
    double rel_Flyr[nlyr];
    read_5c_file ("fpda.atm", &tmp1, &tmp2, &tmp3, &h20, &o3, &nrows);

    //clouds
    double cf = 0.6;
    double *tmp4 = NULL;
    double *tmp5 = NULL; 
    double *tmp6 = NULL;  
    double *tmp7 = NULL;
    double *tmp8 = NULL; 
    double *tmp9 = NULL;  
    double *w0_lw = NULL;
    double *g_lw = NULL;
    
    double *w0_cloud = NULL;
    double *g_cloud = NULL;
    int sw_bands = 0;
    int lw_bands = 0;
    read_5c_file ("rrtm/cldprp/rrtm.sw.int", &tmp4, &tmp5, &tmp6, &w0_cloud, &g_cloud, &sw_bands);
    read_5c_file ("rrtm/cldprp/rrtm.sw.int", &tmp7, &tmp8, &tmp9, &w0_lw, &g_lw, &lw_bands);

    //vertical expansion of cloud
    double tau_ext[sw_bands][nlyr];
    double tau_ext_lw[lw_bands][nlyr];
    
    int cloud_l = 5;
    int cloud_u = 6;
    double tau_ext0 = 13.;
    // tau delta scaling for solar
    for(int insw_bands=0; insw_bands<sw_bands; insw_bands++){
        for(int i=0; i<cloud_l; i++){
	    tau_ext[insw_bands][i] = 0.;
        }
        for(int i=cloud_l; i<=cloud_u; i++){
	    tau_ext[insw_bands][i] = deltatau(tau_ext0,w0_cloud[insw_bands],g_cloud[insw_bands]);
        }
        for(int i=cloud_u+1; i<nlyr; i++){
	    tau_ext[insw_bands][i] = 0.;
        }
    }
    // tau delta sclaing for thermal
    for(int inlw_bands=0; inlw_bands<lw_bands; inlw_bands++){
        for(int i=0; i<cloud_l; i++){
	    tau_ext_lw[inlw_bands][i] = 0.;
        }
        for(int i=cloud_l; i<=cloud_u; i++){
	    tau_ext_lw[inlw_bands][i] = deltatau(tau_ext0,w0_lw[inlw_bands],g_lw[inlw_bands]);
        }
        for(int i=cloud_u+1; i<nlyr; i++){
	    tau_ext_lw[inlw_bands][i] = 0.;
        }
    }
    for(int isw_bands=0; isw_bands>sw_bands; isw_bands++){
	w0_cloud[isw_bands] = deltaw(w0_cloud[isw_bands],g_cloud[isw_bands]);
	g_cloud[isw_bands] = deltag(g_cloud[isw_bands]);
    }
    double g_clear[sw_bands][nlyr];
    double tau_clear[sw_bands][nlyr];
    double w0_clear[sw_bands][nlyr];

    double g_cloudtot[sw_bands][nlyr];
    double tau_cloud[sw_bands][nlyr];
    double w0_cloudtot[sw_bands][nlyr];
    

    //pressure profile
    for (int inlev = 0; inlev < nlev; inlev++) {
        p[inlev] = p0 * (double) inlev / (double) nlyr;
	rel_F[inlev] = h20[inlev]*1e-6*p[inlev]/magnus(tmp3[inlev]);
    }
    for (int inlyr=0; inlyr < nlyr; inlyr++) {
        plyr[inlyr]=(p[inlyr] + p[inlyr+1]) / 2.0;
    }
    double dp = p[1]-p[0];
    for(int inlyr=0; inlyr<nlyr; inlyr++) {
        T[inlyr] = 288. - (nlyr-inlyr-1)*50/(nlyr-1);
	rel_Flyr[inlyr] = (rel_F[inlyr+1] + rel_F[inlyr])/2.;
        o3vmr[inlyr] = 1e-6*(o3[inlyr+1]+o3[inlyr])/2.;
        co2vmr[inlyr] = 400e-6;
        ch4vmr[inlyr] = 1.8e-6;
        n2ovmr[inlyr] = 320e-9;
        o2vmr[inlyr] = .209;
        cfc11vmr[inlyr] = 0;
        cfc12vmr[inlyr] = 0;
        cfc22vmr[inlyr] = 0;
        ccl4vmr[inlyr] = 0;
    }
    
    //time loop	
    while(1==1){
	Ttoa = T[0]; // to compare temperature between 2 time steps
        for (int inlyr=0; inlyr < nlyr; inlyr++) {
                h2ovmr[inlyr] = rel_Flyr[inlyr]*magnus(T[inlyr])/plyr[inlyr];
        }
        //########################## call rrtm every time step ####################################################
	cfpda_rrtm_lw (nlyr, p, T, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr,
		 cfc11vmr,cfc12vmr,cfc22vmr,ccl4vmr,
		 &nbands, &band_lbound,&band_ubound, &wgt_lw, &tau);
	cfpda_rrtm_sw (nlyr,p,T, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr,
		     &nbands_sw, &band_lbound_sw, &band_ubound_sw, &wgt_sw, &dtau_mol_sw, &dtau_ray_sw);
	//file = fopen("kdistribution.asc", "ab");
	/*fprintf(file,"T \n");
	for(int inlyr=0; inlyr<nlyr; inlyr++){
		fprintf(file,"%f  ",T[inlyr]);
	}
	fprintf(file,"\n");
	fprintf(file,"tau \n");*/
	//##################################### cloudless sky ###########################################################
	for (int inbands_sw=0; inbands_sw<nbands_sw; inbands_sw++) {
		for(int inlyr=0; inlyr<nlyr; inlyr++) {
			g_clear[inbands_sw][inlyr] = 0.;
          		tau_clear[inbands_sw][inlyr] = dtau_mol_sw[inbands_sw][inlyr] + dtau_ray_sw[inbands_sw][inlyr];
			w0_clear[inbands_sw][inlyr] = dtau_ray_sw[inbands_sw][inlyr] / tau_clear[inbands_sw][inlyr];
        	}
        }
        kAtmosphere(wgt_lw,band_lbound,band_ubound,nbands,nlyr,nlev,T,tau,Edown,Eup,deltaE_clear);
	solar_rt (deltaE_sw_clear, Ag, nlev, nlyr, wgt_sw, band_lbound_sw, band_ubound_sw, nbands_sw,g_clear,tau_clear,w0_clear);
        //################################################################################################################
	
	//################################### sky with clouds ############################################################
	for (int inbands=0; inbands<nbands; inbands++){
		for (int inlyr=0; inlyr<nlyr; inlyr++){
			tau[inbands][inlyr] = tau[inbands][inlyr]+tau_ext_lw[inbands][inlyr];
		}
	}
	for (int inbands_sw=0; inbands_sw<nbands_sw; inbands_sw++) {
		for(int inlyr=0; inlyr<nlyr; inlyr++) {
			g_cloudtot[inbands_sw][inlyr] = w0_cloud[inbands_sw]*tau_ext[inbands_sw][inlyr]*g_cloud[inbands_sw]/(w0_cloud[inbands_sw]*tau_ext[inbands_sw][inlyr]+dtau_ray_sw[inbands_sw][inlyr]);
          		tau_cloud[inbands_sw][inlyr] = dtau_mol_sw[inbands_sw][inlyr] + dtau_ray_sw[inbands_sw][inlyr] + tau_ext[inbands_sw][inlyr];
			w0_cloudtot[inbands_sw][inlyr] = (w0_cloud[inbands_sw]*tau_ext[inbands_sw][inlyr]+dtau_ray_sw[inbands_sw][inlyr])/tau_cloud[inbands_sw][inlyr];
        	}
        }
        kAtmosphere(wgt_lw,band_lbound,band_ubound,nbands,nlyr,nlev,T,tau,Edown,Eup,deltaE_cloud);
	solar_rt (deltaE_sw_cloud, Ag, nlev, nlyr, wgt_sw, band_lbound_sw, band_ubound_sw, nbands_sw,g_cloudtot,tau_cloud,w0_cloudtot);
	//################################################################################################################

	//################ find dt #######################################################################################
	for(int inlyr=0; inlyr<nlyr; inlyr++){
		deltaT[inlyr] = E2T((1.-cf)*(deltaE_clear[inlyr]+deltaE_sw_clear[inlyr])+cf*(deltaE_cloud[inlyr]+deltaE_sw_cloud[inlyr]),dp);
		if(deltaT[inlyr]<0.){deltaT[inlyr] = -deltaT[inlyr];}
	}
	convection(deltaT,nlyr);
	dt = 1./deltaT[0];
	if(dt > 60*60*12){dt = 60*60*12;}

	//################ new T with calculated dt #####################################################################
	for(int inlyr=0; inlyr<nlyr; inlyr++){
		T[inlyr] += E2T((1.-cf)*(deltaE_clear[inlyr]+deltaE_sw_clear[inlyr])+cf*(deltaE_cloud[inlyr]+deltaE_sw_cloud[inlyr]),dp)*dt;
	}

	//################# convection ##################################################################################
    	for (int inlyr = 0; inlyr < nlyr; inlyr++){
	    theta[inlyr] = T2theta(T[inlyr], plyr[inlyr]);
    	}
    	convection(theta,nlyr);
    	for (int inlyr=0; inlyr < nlyr; inlyr++){
    	    T[inlyr] = theta2T(theta[inlyr],plyr[inlyr]);
    	}
	t += dt;
	//################ convert p to z ##############################################################################
        z[nlev-1] = 0.;
        z[nlev-2] = dp2dz(p[nlev-1]-plyr[nlyr-1],plyr[nlyr-1],T[nlyr-1]);
        for(int i=nlev-3; i>=0; i--){
	    z[i] = z[i+1] + dp2dz(plyr[i+1]-plyr[i],plyr[i],T[i]);
        }
	printf("T(surf) = %6.3f  t(days) = %6.3f \n", T[nlyr-1], t/(60.*60.*24));
	printf("\n");
	//sleep(1);
	//################################### print into file to make plot ###############################################
	/*char buffer[64]; // filename buffer
	snprintf(buffer, sizeof(char) * 64, "./plot/3band_animation/ascii/grey_animation_%i.asc", k);
	file = fopen(buffer, "wb");
	for (int i=0; i < nlyr; i++){
            fprintf(file,"%6.3f %6.3f\n",z[i]*1e-3,T[i]);
        }
	fclose(file );
        */
        //################################################################################################################

	//### if the change in toa temperature is small enough the system is in equilibrium ##############################
	if(((Ttoa-T[0]) < equi && (Ttoa-T[0]) > 0) || ((Ttoa-T[0]) > -equi && (Ttoa-T[0]) < 0)){
	    break;
	}
	
    }
    //ASCII_free_double(tau, nbands);
    //ASCII_free_double(wgt_lw, nbands);
    free(band_lbound);
    free(band_ubound);
    return t;
}
//######################## plots height against temperature for every layer in every time step while running the c program #######################
double plotTimeLoop(double *T, double *theta, double *z, int nlyr, int nlev, int nwvl, double tau_total, double *Edown, double *Eup, double *wvl, double **tau){
    double t = 0.;
    double p0 = 1000.0; //hPa
    double equi = 1e-5;

    double p[nlev],plyr[nlyr];
    double Tnlev[nlev],deltaT[nlyr],dt,Tsurf;
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
    FILE *fp;    

    int count = 0;
    for (double tau_total=tau_min; tau_total<=tau_max; tau_total += tau_step){
        t = timeLoop(T,theta,z,nlyr,nlev,nwvl,tau_total,Edown,Eup,wvl,taulbl);
        tau[count] = tau_total;
        Tequ[count] = T[nlyr-1];
	count++;
    }
    fp = fopen("./plot/tau/3band.asc", "w+");
    for (int i=0; i < count; i++){
        fprintf(fp,"%6.3f %6.3f\n",Tequ[i],tau[i]);
    }
    return t;
}
int main()
{ 
    //####################### for grey and 3-band atmospheres ############################################################
    /*int nlyr = 30;
    int nlev = nlyr + 1;
    int nwvl = 3;
    double z[nlev];
    double T[nlyr],theta[nlyr];
    double Eup[nlev],Edown[nlev];*/
    //####################################################################################################################
  
    //##################### for plotting different taus in grey and 3-band atmospheres ###################################
    /*double tau_min = 0.;
    double tau_max = 2.0;
    int tau_size = (tau_max-tau_min)/0.1 + 1;*/
    //####################################################################################################################
   
    //################################## for line by line atmosphere #####################################################
    /*double t;
    int nwvlco2 = 0;
    int tmp1 = 0;
    int nyco2 = 0;
    int tmp2 = 0;
    double tau_total = 1.9;
    double *wvnco2 = NULL; 
    double **tauco2 = NULL; 
    double *tmp3 = NULL; 
    double **tauh2o = NULL;
    double **tauch4 = NULL;
    double **taun2o = NULL;
    double **tauo3 = NULL;
   
    ASCII_file2xy2D ("lbl.arts/lbl.co2.asc", &nwvlco2, &nyco2, &wvnco2, &tauco2);
    ASCII_file2xy2D ("lbl.arts/lbl.h2o.asc", &tmp1, &tmp2, &tmp3, &tauh2o);
    ASCII_file2xy2D ("lbl.arts/lbl.ch4.asc", &tmp1, &tmp2, &tmp3, &tauch4);
    ASCII_file2xy2D ("lbl.arts/lbl.n2o.asc", &tmp1, &tmp2, &tmp3, &taun2o);
    ASCII_file2xy2D ("lbl.arts/lbl.o3.asc", &tmp1, &tmp2, &tmp3, &tauo3);
   
    int nlyrlbl = nyco2;
    int nlevlbl = nlyrlbl + 1;
    double zlbl[nlevlbl];
    double Tlbl[nlyrlbl],thetalbl[nlyrlbl];
    double Euplbl[nlevlbl],Edownlbl[nlevlbl];
    //add up tau's from h2o and co2
    for (int inwvl = 0; inwvl<nwvlco2; inwvl++){
	for (int inlyr = 0; inlyr<nlyrlbl; inlyr++){
		tauco2[inwvl][inlyr] = tauco2[inwvl][inlyr] + tauh2o[inwvl][inlyr] + tauch4[inwvl][inlyr] + taun2o[inwvl][inlyr] + tauo3[inwvl][inlyr];
	}
    }*/
    //###################################################################################################################
   
    //################################## for k-distribution atmosphere ##################################################
    double t;
    int knlyr = 20;
    int knlev = knlyr + 1;
    //###################################################################################################################

    //############################ call different functions #############################################################
    t = kTimeLoop(knlyr,knlev); //time loop for k distribution
    //t = timeLoop(T,theta,z,nlyr,nlev,nwvl,tau_total,Edown,Eup,wvnco2,tauco2); //for grey and window atmosphere
    //t = timeLoop(Tlbl,thetalbl,zlbl,nlyrlbl,nlevlbl,nwvlco2,tau_total,Edownlbl,Euplbl,wvnco2,tauco2); //for line by line atmosphere
    //t = plotTimeLoop(T,theta,z,nlyr,nlev,nwvl,tau_total,Edown,Eup,wvnco2,tauco2);
    //t = plotTimeLoop(Tlbl,thetalbl,zlbl,nlyrlbl,nlevlbl,nwvlco2,tau_total,Edownlbl,Euplbl,wvnco2,tauco2);
    //t = plotT2tau(T,theta,z,nlyr,nlev,nwvl,tau_min,tau_max,tau_size,Edown,Eup,wvnco2,tauco2);
    //###################################################################################################################
    
    //################################### print the output into a file ##################################################
    /*FILE *fp;
    fp = fopen("profile_lbl.txt", "ab");
    //fprintf(fp, "after %6.2f years\n",t/(365.*24*60*60)); //time until equilibrium
    //fprintf(fp,"Eup[TOA] = %6.3f\n",Eup[0]);             //Eup at top of atmosphere, should be Eearth when equilibrium is reached
    //fprintf(fp,"z[km], T[K]:\n");                        //height and temperature profile
    for (int i=0; i < nlyr; i++){
        fprintf(fp,"%6.3f %6.3f \n",zlbl[i]*1e-3,Tlbl[i]);
    }
    fclose(fp);*/
    //###################################################################################################################

    /*ASCII_free_double(tauh2o, nwvlco2);
    ASCII_free_double(tauco2, nwvlh2o);
    free(wvnco2);
    free(wvnh2o);*/
}
