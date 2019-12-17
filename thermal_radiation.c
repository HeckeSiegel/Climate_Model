#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define pi 3.141592
#define h 6.62607015e-34 //Js
#define c 2.99792458e8 // m/s
#define kB 1.38064852e-23 // m^2kgs^-2K^-1
#define sigma 5.670374e-8
#define Eearth 235.

double B_gray(double T){
    return sigma*pow(T,4)/pi;
}
double B(double wvl, double T){
    return 2*h*pow(c,2)/pow(wvl,5) * 1./(exp(h*c/(wvl*kB*T))-1.);
}
double B_int(double lambda_1, double lambda_2, double T){
    if(lambda_2-lambda_1 < 1e-6){
	return B(lambda_1+(lambda_2-lambda_1)/2.,T)*(lambda_2-lambda_1);
    }
    else{
    	int N = (int)ceil((lambda_2-lambda_1))/1e-6;
	double dwvl;
	if(N>50){
		N = 50;
		dwvl = (lambda_2-lambda_1)/N;
	}
	else{
    		dwvl = (lambda_2-lambda_1)/N;
	}
    	double wvl = lambda_1+dwvl/2.;
    	double res = 0.;
    	//integration
    	for(int i=0; i<N; i++){
        	res += B(wvl+dwvl*i,T)*dwvl;
    	}
    	return res;
    }
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
    deltaE[nlyr-1] = Edown[nlyr-1] - Eup[nlyr-1];
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
    for (int iwvl=0; iwvl<nwvl; iwvl++){
        B_surface = B_int(wvlband[iwvl][0],wvlband[iwvl][1],T[nlyr-1]);
        for (int inlyr = 0; inlyr< nlyr; inlyr++){
            B_layer[inlyr] = B_int(wvlband[iwvl][0],wvlband[iwvl][1],T[inlyr]);
	    tmp_tau[inlyr] = tau[iwvl][inlyr];
        }
    	schwarzschild(nlev,tmp_tau,B_layer,B_surface,tmp_Edown,tmp_Eup);
	for(int inlev=0; inlev<nlev; inlev++){
	    Eup[inlev] += tmp_Eup[inlev];
	    Edown[inlev] += tmp_Edown[inlev];
        }
    }
    dE(deltaE,Edown,Eup,nlyr);
}
void multiBandAtmosphere(double *wvl,int nwvl, int nlyr, int nlev, double *T, double **tau, double *Edown, double *Eup, double *deltaE){
    double tmp_tau[nlyr];
    double dwvl,central_wvl, B_surface, B_layer[nlyr];
    double tmp_Edown[nlev], tmp_Eup[nlev];

    for (int inlev=0; inlev<nlev; inlev++){
	Eup[inlev] = 0.;
	Edown[inlev] = 0.;
    }
    for (int iwvl=0; iwvl<(nwvl-1); iwvl++){
	dwvl = (wvl[iwvl+1]-wvl[iwvl])*1e-9;
	central_wvl = wvl[iwvl]*1e-9 + dwvl/2.;
        B_surface = B(central_wvl,T[nlyr-1])*dwvl;
        for (int inlyr = 0; inlyr< nlyr; inlyr++){
            B_layer[inlyr] = B(central_wvl,T[inlyr])*dwvl;
	    tmp_tau[inlyr] = tau[iwvl][inlyr];
        }
    	schwarzschild(nlev,tmp_tau,B_layer,B_surface,tmp_Edown,tmp_Eup);
	for(int inlev=0; inlev<nlev; inlev++){
	    Eup[inlev] += tmp_Eup[inlev];
	    Edown[inlev] += tmp_Edown[inlev];
        }
    }
    dE(deltaE,Edown,Eup,nlyr);
}
void kAtmosphere(double **wgt_lw, double *band_lbound, double *band_ubound, int nbands, int nlyr, int nlev, double *T, double **tau, double *Edown, double *Eup, double *deltaE){
    double tmp_tau[nlyr];
    double B_surface, B_layer[nlyr], tmp_B_surface_int, tmp_B_layer_int[nlyr];
    double tmp_Edown[nlev], tmp_Eup[nlev];

    for (int inlev=0; inlev<nlev; inlev++){
	Eup[inlev] = 0.;
	Edown[inlev] = 0.;
    }

    tmp_B_surface_int = B_int(1e-2/band_ubound[0],1e-2/band_lbound[0],T[nlyr-1]);
    B_surface = tmp_B_surface_int*wgt_lw[0][nlyr-1]; 
    for (int inlyr = 0; inlyr< nlyr; inlyr++){
	tmp_B_layer_int[inlyr] = B_int(1e-2/band_ubound[0],1e-2/band_lbound[0],T[inlyr]);
        B_layer[inlyr] = tmp_B_layer_int[inlyr]*wgt_lw[0][inlyr];
	tmp_tau[inlyr] = tau[0][inlyr];
    }
    schwarzschild(nlev,tmp_tau,B_layer,B_surface,tmp_Edown,tmp_Eup);
    for(int inlev=0; inlev<nlev; inlev++){
	Eup[inlev] += tmp_Eup[inlev];
	Edown[inlev] += tmp_Edown[inlev];
    }

    for (int inbands=1; inbands<nbands; inbands++){
	if(band_ubound[inbands-1]==band_ubound[inbands] && band_lbound[inbands-1]==band_lbound[inbands])
	{
	    B_surface = tmp_B_surface_int*wgt_lw[inbands][nlyr-1];
	    for (int inlyr = 0; inlyr< nlyr; inlyr++){
            B_layer[inlyr] = tmp_B_layer_int[inlyr]*wgt_lw[inbands][inlyr];
	    tmp_tau[inlyr] = tau[inbands][inlyr];
            }
	}
	else{
	    tmp_B_surface_int = B_int(1e-2/band_ubound[inbands],1e-2/band_lbound[inbands],T[nlyr-1]);
            B_surface = tmp_B_surface_int*wgt_lw[inbands][nlyr-1];
            for (int inlyr = 0; inlyr< nlyr; inlyr++){
	    	tmp_B_layer_int[inlyr] = B_int(1e-2/band_ubound[inbands],1e-2/band_lbound[inbands],T[inlyr]);
            	B_layer[inlyr] = tmp_B_layer_int[inlyr]*wgt_lw[inbands][inlyr];
	    	tmp_tau[inlyr] = tau[inbands][inlyr];
            }
	}
    	schwarzschild(nlev,tmp_tau,B_layer,B_surface,tmp_Edown,tmp_Eup);
	for(int inlev=0; inlev<nlev; inlev++){
	    Eup[inlev] += tmp_Eup[inlev];
	    Edown[inlev] += tmp_Edown[inlev];
        }
    }
    dE(deltaE,Edown,Eup,nlyr);
}
/* calculate layer properties t, r, rdir, sdir, and tdir from            */
/* layer optical thickness dtau, asymmetry parameter g,                  */
/* single scattering albedo omega0, and cosine of solar zenith angle mu0 */

void eddington_v2 (double dtau, double g, double omega0, double mu0,
		   double *t, double *r, double *rdir, double *sdir, double *tdir)
{
  double alpha1=0, alpha2=0, alpha3=0, alpha4=0, alpha5=0, alpha6=0;
  double a11=0, a12=0, a13=0, a23=0, a33=0;
  double lambda=0, b=0, A=0;
  double denom=0;

  /* first, avoid omega0=1 because of instability */
  if (omega0 > 0.999999)
    omega0=0.999999;
  if (dtau>100.0){
    dtau = 100.0;
  }
  alpha1= (1.0-omega0)+0.75*(1.0-omega0*g);
  alpha2=-(1.0-omega0)+0.75*(1.0-omega0*g);
  
  lambda=sqrt(alpha1*alpha1-alpha2*alpha2);
  
  A=1.0/(alpha2/(alpha1-lambda)*exp(lambda*dtau)-alpha2/(alpha1+lambda)*exp(-lambda*dtau));
  
  a11=A*2.0*lambda/alpha2;
  a12=A*(exp(lambda*dtau)-exp(-lambda*dtau));
  
  b=0.5-0.75*g*mu0;
  alpha3=-omega0*b; 
  alpha4=omega0*(1-b);
  
  denom = (1.0/mu0/mu0-lambda*lambda);
  alpha5=((alpha1-1.0/mu0)*alpha3-alpha2*alpha4)/denom;
  alpha6=(alpha2*alpha3-(alpha1+1.0/mu0)*alpha4)/denom;
  
  a33=exp(-dtau/mu0);
  
  a13=alpha5*(1.0-(a11)*(a33))-alpha6*(a12);
  a23=-(a12)*alpha5*(a33)+alpha6*((a33)-(a11));

  *t    = a11;
  *r    = a12;
  *tdir = a33;
  *rdir = a13 / mu0;
  *sdir = a23 / mu0;
}
void solar_rt (double *deltaE, double Ag, int nlev, int nlyr, double **tau, double g, double **omega0, double *wgt_sw, double *band_lbound, double *band_ubound, int nbands){
	double r[nlyr], t[nlyr], tdir[nlyr], rdir[nlyr], sdir[nlyr]; 
	double R[nlyr], T[nlyr], Tdir[nlyr], Sdir[nlyr];
	double Eet = 1361.6; // W/m^2
	double mu = 0.25; // cos(60°)/2
	double Edown[nlev], Eup[nlev], Edir[nlev], E_nsdown[nlev], E_nsup[nlev];
	double Edown_total[nlev], Eup_total[nlev];

	for (int inlev=0; inlev<nlev; inlev++){
		Eup_total[inlev] = 0.;
		Edown_total[inlev] = 0.;
    	}
	//calculate deltaE's for all wavelength bands
	for (int inbands=0; inbands<nbands; inbands++){
	for (int inlyr=0; inlyr<nlyr; inlyr++){
		eddington_v2 (tau[inbands][inlyr],g,omega0[inbands][inlyr],mu,&(t[inlyr]),&(r[inlyr]),&(rdir[inlyr]),&(sdir[inlyr]),&(tdir[inlyr]));
	}
	// 1. R,Tdir,Sdir
	R[0] = r[0];
	T[0] = t[0];
	Tdir[0] = tdir[0];
	Sdir[0] =  (t[1]*sdir[0]+tdir[0]*rdir[1]*r[0]*t[1])/(1. - r[0]*r[1]);
	for (int inlyr=1; inlyr<nlyr; inlyr++){
		R[inlyr] = r[inlyr] + R[inlyr-1]*t[inlyr]*t[inlyr]/(1. - R[inlyr-1]*r[inlyr]);
		T[inlyr] = T[inlyr-1]*t[inlyr]/(1. - R[inlyr-1]*r[inlyr]);
		Tdir[inlyr] *= tdir[inlyr];
		Sdir[inlyr] = Tdir[inlyr-1]*sdir[inlyr] + (t[inlyr]*Sdir[inlyr-1] + Tdir[inlyr-1]*rdir[inlyr]*R[inlyr-1]*t[inlyr])/(1. - R[inlyr-1]*r[inlyr]);
	}
	// 2. irradiance at surface
	// irradiance without source term (ns = no source)
	E_nsdown[0] = Eet*mu;
	E_nsdown[nlev-1] = T[nlev-2]*E_nsdown[0]/(1. - R[nlev-2]*Ag);
	E_nsup[nlev-1] = Ag*E_nsdown[nlev-1];
	for (int inlev=nlev-2; inlev>0; inlev--){
		E_nsdown[inlev] = (T[inlev-1]*E_nsdown[0] + R[inlev-1]*t[inlev-1]*E_nsup[inlev+1])/(1. - R[inlev-1]*r[inlev]);
		E_nsup[inlev] = (t[inlev]*E_nsup[inlev+1] + T[inlev-1]*r[inlev]*E_nsdown[0])/(1. - R[inlev-1]*r[inlev]);
	} 
	E_nsup[0] = t[0]*E_nsup[1] + r[0]*E_nsdown[0];
	// direct irradiance
	Edir[0] = Eet*mu*wgt_sw[inbands];
	for (int inlev=1; inlev<nlev; inlev++){
		Edir[inlev] = Edir[inlev-1]*tdir[inlev-1];
	}
	// total irradiance at surface
	Edown[nlev-1] = (Sdir[nlev-1] + Tdir[nlev-2]*R[nlev-2]*Ag)*Edir[0]/(1. - R[nlev-2]*Ag);
	Eup[nlev-1] = Ag*(Edown[nlev-1] + Tdir[nlev-2]*Edir[0]);
	
	// 3. irradiances at all other levels
	for (int inlev=nlev-2; inlev>0; inlev--){
		Edown[inlev] = E_nsdown[inlev] + (Edir[0]*Sdir[inlev] + Edir[inlev]*rdir[inlev]*R[inlev-1])/(1. - R[inlev-1]*r[inlev]);
		Eup[inlev] = E_nsup[inlev] + (Edir[0]*Sdir[inlev]*rdir[inlev] + Edir[inlev]*rdir[inlev])/(1. - R[inlev-1]*r[inlev]);
	}
	Edown[0] = 0.;
	Eup[0] = t[0]*Eup[1] + rdir[0]*Edir[0];
	for (int inlev=0; inlev<nlev; inlev++){
		Eup_total[inlev] += Eup[inlev];
		Edown_total[inlev] += Edown[inlev];
    	}
	}
	
	dE(deltaE,Edown_total,Eup_total,nlyr);
}
