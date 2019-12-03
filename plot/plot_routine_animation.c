#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include "gnuplot_i.h"
#include "ascii.h"

int main(){

    	int ny = 0;
    	int nx = 0;
    	double *z = NULL; 
    	double **T0 = NULL; 
	

	for(int k = 0; k<1184; k=k+10){

	char buffer[64]; // filename buffer
	snprintf(buffer, sizeof(char) * 64, "./3band_animation/ascii/grey_animation_%i.asc", k);

    	ASCII_file2xy2D (buffer, &nx, &ny, &z, &T0);
	double T[nx];

	for(int inx = 0; inx<nx; inx++){
		T[inx] = T0[inx][0];
	}	
	
	gnuplot_ctrl *g1;
    	g1 = gnuplot_init();	

	//gnuplot_resetplot  (g1);  /* start with new plot rather than plotting into exisiting one */
    	gnuplot_setstyle   (g1, "lines");      /* draw lines and points */
    	gnuplot_set_xlabel (g1, "temperature [K]");  /* xaxis label */
    	gnuplot_set_ylabel (g1, "height [km]");    /* yaxis label */

	char buffer2[64];
	gnuplot_cmd(g1,"set terminal png" );
	snprintf(buffer2, sizeof(char) * 64, "set output './3band_animation/png/grey_animation_%i.png'", k);
	gnuplot_cmd(g1, buffer2);
    	/* plot temperature T as function of z and label with temperature */
    	gnuplot_plot_xy   (g1, T, z, nx, "Temperature") ;

    	sleep(1);	 
    	/* close plot */
    	gnuplot_close (g1) ;

	
	}
        ASCII_free_double(T0, nx);
        free(z);
	return 0;
}
