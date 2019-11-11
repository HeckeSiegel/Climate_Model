#include <stdio.h>
#include <unistd.h>

/* gnuplot_i libary */
#include "gnuplot_i.h"


int main (int argc, char **argv)
{
  /* temperature profile consists of only four points */
  int nlyr=4;
  double z[4] = {0, 10, 15, 50};
  double T[4] = {288, 223, 223, 240};  

  int tcounter=0;

  gnuplot_ctrl *g1;
  g1 = gnuplot_init();
  
  /* infinite loop over temperature profile */
  while (1==1) {
    /* number of time steps */
    tcounter++;

    T[0] = T[0] - 0.1;  /* decrease surface temperature by 0.1K */
    printf ("T[0] = %f\n", T[0]);

    /* plot every 10th time step */
    if (tcounter%10 == 0) {
      
      gnuplot_resetplot  (g1);  /* start with new plot rather than plotting into exisiting one */
      gnuplot_setstyle   (g1, "linespoints");      /* draw lines and points */
      gnuplot_set_xlabel (g1, "temperature [K]");  /* xaxis label */
      gnuplot_set_ylabel (g1, "altitude [km]");    /* yaxis label */
      
      /* plot temperature T as function of z and label with temperature */
      gnuplot_plot_xy   (g1, T, z, nlyr, "Temperature") ;
      sleep(1); /* wait a second */
    }
  }

  /* close plot */
  gnuplot_close (g1) ;

  return 0;
}
