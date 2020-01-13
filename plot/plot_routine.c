#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include "gnuplot_i.h"
#include "ascii.h"

int main(){


	/*int nxzero = 0;
	int nxnormal = 0;
	int nxdouble = 0;
	int nxall = 0;*/
	int nxlbl = 0;
	int tmpy = 0;

	/*double *Tzero = NULL;
	double *Tnormal = NULL;
	double *Tdouble = NULL;
	double *Tall = NULL;*/
	double *zlbl = NULL;

	/*double **tzero0 = NULL; 
	double **tnormal0 = NULL;
	double **tdouble0 = NULL;
	double **tall0 = NULL;*/
	double **Tlbl0 = NULL;
	
	

	/*ASCII_file2xy2D ("lbl_zero_co2_new.asc", &nxzero, &tmpy, &Tzero, &tzero0);
	ASCII_file2xy2D ("lbl_normal_co2_new.asc", &nxnormal, &tmpy, &Tnormal, &tnormal0);
	ASCII_file2xy2D ("lbl_double_co2_new.asc", &nxdouble, &tmpy, &Tdouble, &tdouble0);
	ASCII_file2xy2D ("lbl_zero_all_new.asc", &nxall, &tmpy, &Tall, &tall0);*/
	ASCII_file2xy2D ("profile_lbl.txt", &nxlbl, &tmpy, &zlbl, &Tlbl0);

	/*double tzero[nxzero];
	double tnormal[nxnormal];
	double tdouble[nxdouble];
	double tall[nxall];*/
	double Tlbl[nxlbl];

	/*for(int inx1 = 0; inx1<nxzero; inx1++){
		tzero[inx1] = tzero0[inx1][0];
	}	
	for(int inx2 = 0; inx2<nxnormal; inx2++){
		tnormal[inx2] = tnormal0[inx2][0];
	}
	for(int inx3 = 0; inx3<nxdouble; inx3++){
		tdouble[inx3] = tdouble0[inx3][0];
	}
	for(int inx4 = 0; inx4<nxall; inx4++){
		tall[inx4] = tall0[inx4][0];
	}*/
	for(int inx5 = 0; inx5<nxlbl; inx5++){
		Tlbl[inx5] = Tlbl0[inx5][0];
	}

	gnuplot_ctrl *g1;
    	g1 = gnuplot_init();	

	//gnuplot_resetplot  (g1);  /* start with new plot rather than plotting into exisiting one */
    	gnuplot_setstyle   (g1, "lines");      /* draw lines and points */
    	gnuplot_set_xlabel (g1, "height [km]");  /* xaxis label */
    	gnuplot_set_ylabel (g1, "temperature [K]");    /* yaxis label */
    	
    	/*gnuplot_plot_xy   (g1, tall, Tall, nxall, "no greenhouse gases") ;
	gnuplot_plot_xy   (g1, tzero, Tzero, nxzero, "CO2 = 0 ppm") ;
	gnuplot_plot_xy   (g1, tnormal, Tnormal, nxnormal, "CO2 = 400 ppm") ;
	gnuplot_plot_xy   (g1, tdouble, Tdouble, nxdouble, "CO2 = 800 ppm") ;*/
	gnuplot_plot_xy   (g1, Tlbl, zlbl, nxlbl," ");

    	sleep(50);	 
    	/* close plot */
    	//gnuplot_close (g1) ;


	/*ASCII_free_double(tall0, nxall);
	ASCII_free_double(tnormal0, nxnormal);
	ASCII_free_double(tdouble0, nxdouble);
	ASCII_free_double(tzero0, nxzero);*/
	ASCII_free_double(Tlbl0, nxlbl);
	free(zlbl);
        /*free(Tall);
	free(Tnormal);
	free(Tdouble);
	free(Tzero);*/
	return 0;
}
