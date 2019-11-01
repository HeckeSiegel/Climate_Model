#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define pi 3.14

double alpha(double mu, double tau){
    return 1. - exp(-tau/mu);
}
void schwarzschild(int nlev, double *tau, double *B_layer, double B_surface, double *Edown, double *Eup){
    int Nmu = 10;
    double dmu = 1./Nmu;
    double mu[Nmu],Ldown0[Nmu], Ldown1[Nmu];
    //Edown
    for (int iNmu = 0; iNmu < Nmu; iNmu++){
        Ldown0[iNmu] = 0.;
	mu[iNmu] = 0.05 + iNmu*dmu;
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
    for (int i = 0; i<nlyr; i++){
	deltaE[i] = Edown[i] + Eup[i+1] - Eup[i] - Edown[i+1];
	}
}
