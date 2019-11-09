#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define pi 3.141592
#define h 6.62607015e-34 //Js
#define c 2.99792458e8 // m/s
#define kB 1.38064852e-23 // m^2kgs^-2K^-1
#define sigma 5.67e-8
#define Eearth 240.4
#define Evenus 261.

double B_gray(double T){
    return sigma*pow(T,4)/pi;
}
double B(double wvl, double T){
    return 2*h*pow(c,2)/pow(wvl,5) * 1./(exp(h*c/(wvl*kB*T))-1.);
}
double B_int(double lambda_1, double lambda_2, double T){
    double dwvl = 1e-7;
    int N = (lambda_2-lambda_1)/dwvl;
    double wvl = lambda_1+dwvl/2.;
    double res = 0.;
    //integration
    for(int i=0; i<N; i++){
        res += B(wvl+dwvl*i,T)*dwvl;
    }
    return res;
}
double alpha(double mu, double tau){
    return 1. - exp(-tau/mu);
}
void schwarzschild(int nlev, double *tau, double *B_layer, double B_surface, double *Edown, double *Eup){
    int Nmu = 20;
    double dmu = 1./Nmu;
    double mu[Nmu],Ldown0[Nmu], Ldown1[Nmu];
    //Edown
    for (int iNmu = 0; iNmu < Nmu; iNmu++){
        Ldown0[iNmu] = 0.;
	mu[iNmu] = dmu/2. + iNmu*dmu;
    }
    Edown[0] = 0.;
    for (int inlev = 1; inlev < nlev; inlev++){
        double tmp = 0;
        for (int iNmu = 0; iNmu < Nmu; iNmu++){
	    Ldown1[iNmu] = (1 - alpha(mu[iNmu],tau[inlev-1]))*Ldown0[iNmu] + alpha(mu[iNmu],tau[inlev-1])*B_layer[inlev-1];
            tmp += Ldown1[iNmu]*mu[iNmu];
            Ldown0[iNmu] = Ldown1[iNmu];
        }
    Edown[inlev] = tmp*2.*pi*dmu;
    }
    //Eup
    double Lup0[Nmu], Lup1[Nmu];
    for (int iNmu = 0; iNmu<Nmu; iNmu++){
        Lup0[iNmu] = B_surface;
    }
    Eup[nlev-1] = B_surface*pi;
    for (int inlev = nlev-2; inlev >= 0; inlev--){
        double tmp = 0;
        for (int iNmu = 0; iNmu < Nmu; iNmu++){
	    Lup1[iNmu] = (1 - alpha(mu[iNmu],tau[inlev]))*Lup0[iNmu] + alpha(mu[iNmu],tau[inlev])*B_layer[inlev];
            tmp += Lup1[iNmu]*dmu*mu[iNmu];
	    Lup0[iNmu] = Lup1[iNmu];
        }
    Eup[inlev] = tmp*2.*pi;
    }
}
void dE(double *deltaE, double *Edown, double *Eup, int nlyr){
    for (int i = 0; i<nlyr-1; i++){
	deltaE[i] = Edown[i] + Eup[i+1] - Eup[i] - Edown[i+1];
	}
    deltaE[nlyr-1] = Edown[nlyr-1] - Eup[nlyr-1] + Eearth;
}
void oneBandAtmosphere(double *T, int nlev, int nlyr, double tau_total,double *Edown, double *Eup, double *deltaE){
    double tau[nlyr];
    double B_surface, B_layer[nlyr];
    double dtau = tau_total/nlyr;
    for (int inlyr = 0; inlyr< nlyr; inlyr++){
        tau[inlyr] = dtau;
    }
    B_surface = B_gray(T[nlyr-1]);
    for(int inlyr=0; inlyr<nlyr; inlyr++){
	B_layer[inlyr] = B_gray(T[inlyr]);
    }
    schwarzschild(nlev,tau,B_layer,B_surface,Edown,Eup);
    dE(deltaE,Edown,Eup,nlyr);
}
void threeBandAtmosphere(int nwvl, int nlyr, int nlev, double *T, double tau_total, double *Edown, double *Eup, double *deltaE){
    double tau[nwvl][nlyr], tmp_tau[nlyr];
    double B_surface, B_layer[nlyr];
    double dtau = tau_total/nlyr;
    double tmp_Edown[nlev], tmp_Eup[nlev];
    //double tmp_B_surface,tmp_B_layerF,tmp_B_layerN;

    for (int inlyr = 0; inlyr< nlyr; inlyr++){
	for (int inwvl = 0; inwvl<nwvl; inwvl++){
        tau[inwvl][inlyr] = dtau;
	}
    }
    for (int inlyr = 0; inlyr< nlyr; inlyr++){
	tau[1][inlyr] = 0.;
    }
    double wvlband[3][2] = {{1e-6,8e-6},{8e-6,12e-6},{12e-6,1e-3}};
    for (int inlev=0; inlev<nlev; inlev++){
	Eup[inlev] = 0.;
	Edown[inlev] = 0.;
    }
    /*tmp_B_surface = 0.;
    tmp_B_layerF = 0.;
    tmp_B_layerN = 0.;*/
    for (int iwvl=0; iwvl<nwvl; iwvl++){
        B_surface = B_int(wvlband[iwvl][0],wvlband[iwvl][1],T[nlyr-1]);
        for (int inlyr = 0; inlyr< nlyr; inlyr++){
            B_layer[inlyr] = B_int(wvlband[iwvl][0],wvlband[iwvl][1],T[inlyr]);
	    tmp_tau[inlyr] = tau[iwvl][inlyr];
        }
	/*tmp_B_surface += B_surface;
	tmp_B_layerF += B_layer[5];
	tmp_B_layerN += B_layer[0];*/
    	schwarzschild(nlev,tmp_tau,B_layer,B_surface,tmp_Edown,tmp_Eup);
	for(int inlev=0; inlev<nlev; inlev++){
	    Eup[inlev] += tmp_Eup[inlev];
	    Edown[inlev] += tmp_Edown[inlev];
        }
    }
    dE(deltaE,Edown,Eup,nlyr);
    //printf("%6.3f %6.3f %6.3f\n",tmp_B_surface,tmp_B_layerF,tmp_B_layerN);
}
