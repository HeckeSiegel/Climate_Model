#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define pi 3.14

double alpha(double mu, double tau){
    return 1. - exp(-tau/mu);
}
void schwarzschild(int nlev, double *tau, double *B_layer, double B_surface, double *Edown, double *Eup){
    int Nmu = 10;
    //mu profile
    double dmu = 1./Nmu;
    double mu[Nmu];
    for (int i = 0; i < Nmu; i++){
        mu[i] = 0.05 + i*dmu;
    }
    //Edown
    double Ldown[nlev][Nmu];
    for (int iN = 0; iN<Nmu; iN++){
        Ldown[0][iN] = 0.;
    }
    for (int inlev = 1; inlev < nlev; inlev++){
        for (int iN = 0; iN < Nmu; iN++){
            Ldown[inlev][iN] = (1 - alpha(mu[iN],tau[inlev-1]))*Ldown[inlev-1][iN] + alpha(mu[iN],tau[inlev-1])*B_layer[inlev-1];
        }
    }
    for (int inlev = 0; inlev < nlev; inlev++){
        double tmp = 0;
        for (int iN = 0; iN < Nmu; iN++){
            tmp += Ldown[inlev][iN]*mu[iN];
        }
    Edown[inlev] = tmp*2.*pi*dmu;
    }
    //Eup
    double Lup[nlev][Nmu];
    for (int iN = 0; iN<Nmu; iN++){
        Lup[nlev-1][iN] = B_surface;
    }
    for (int inlev = nlev-2; inlev >= 0; inlev--){
        for (int iN = 0; iN < Nmu; iN++){
            Lup[inlev][iN] = (1 - alpha(mu[iN],tau[inlev]))*Lup[inlev+1][iN] + alpha(mu[iN],tau[inlev])*B_layer[inlev];
        }
    }
    for (int inlev = 0; inlev < nlev; inlev++){
        double tmp = 0;
        for (int iN = 0; iN < Nmu; iN++){
            tmp += Lup[inlev][iN]*dmu*mu[iN];
        }
    Eup[inlev] = tmp*2.*pi;
    }
}
void dE(double *deltaE, double *Edown, double *Eup, int nlyr){
    for (int i = 0; i<nlyr; i++){
	deltaE[i] = Edown[i] + Eup[i+1] - Eup[i] - Edown[i+1];
	}
}
