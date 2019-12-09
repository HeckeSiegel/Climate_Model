#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include "gnuplot_i.h"
#include "ascii.h"

int main(){

    	int ny = 0;
    	int nx = 0;
    	double *T288 = NULL; 
	double *Tgrey = NULL; 
	double *T3band = NULL; 
    	double **tau0 = NULL; 
	
    	ASCII_file2xy2D ("./tau/288.asc", &nx, &ny, &T288, &tau0);
	ASCII_file2xy2D ("./tau/grey.asc", &nx, &ny, &Tgrey, &tau0);
	ASCII_file2xy2D ("./tau/3band.asc", &nx, &ny, &T3band, &tau0);
	double tau[nx];

	for(int inx = 0; inx<nx; inx++){
		tau[inx] = tau0[inx][0];
	}	

	gnuplot_ctrl *g1;
    	g1 = gnuplot_init();	

	//gnuplot_resetplot  (g1);  /* start with new plot rather than plotting into exisiting one */
    	gnuplot_setstyle   (g1, "lines");      /* draw lines and points */
    	gnuplot_set_xlabel (g1, "total optical thickness of the atmosphere");  /* xaxis label */
    	gnuplot_set_ylabel (g1, "equilibrium surface temperature [K]");    /* yaxis label */
    	
    	gnuplot_plot_xy   (g1, tau, Tgrey, nx, "grey atmosphere") ;
	gnuplot_plot_xy   (g1, tau, T3band, nx, "3 band atmosphere") ;
	gnuplot_plot_xy   (g1, tau, T288, nx, "average surface temperature earth") ;

    	sleep(50);	 
    	/* close plot */
    	//gnuplot_close (g1) ;


	ASCII_free_double(tau0, nx);
        free(T288);
	free(Tgrey);
	free(T3band);
	return 0;
}
