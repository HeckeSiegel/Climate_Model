#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include "gnuplot_i.h"
#include "ascii.h"

int main(){

    	int ny = 0;
    	int nx = 0;
    	double *T = NULL; 
    	double **t0 = NULL; 
	
    	ASCII_file2xy2D ("lbl_normal_co2.asc", &nx, &ny, &T, &t0);
	double t[nx];

	for(int inx = 0; inx<nx; inx++){
		t[inx] = t0[inx][0];
	}	

	gnuplot_ctrl *g1;
    	g1 = gnuplot_init();	

	//gnuplot_resetplot  (g1);  /* start with new plot rather than plotting into exisiting one */
    	gnuplot_setstyle   (g1, "lines");      /* draw lines and points */
    	gnuplot_set_xlabel (g1, "time (days)");  /* xaxis label */
    	gnuplot_set_ylabel (g1, "surface temperature [K]");    /* yaxis label */
    	/* plot temperature T as function of z and label with temperature */
    	gnuplot_plot_xy   (g1, t, T, nx, "Temperature") ;
    
    	sleep(50);	 
    	/* close plot */
    	//gnuplot_close (g1) ;


	ASCII_free_double(t0, nx);
        free(T);
	return 0;
}
