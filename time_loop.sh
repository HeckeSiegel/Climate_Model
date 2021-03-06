CC=gcc

CFLAGS='-O3 -ffast-math -fstrict-aliasing -g -Wall' # -fsanitize=address

DLFLAGS='-lm -lfpda_rrtm_lw -lfpda_rrtm_sw'

$CC -c $CFLAGS thermal_radiation.c

$CC -c $CFLAGS -I./gnuplot/gnuplot_i-2.10/src/ -I./lbl.arts time_loop.c -I./rrtm/install/include/ 

gfortran $CFLAGS -L./rrtm/install/lib/ -o time_loop time_loop.o thermal_radiation.o gnuplot/gnuplot_i-2.10/gnuplot_i.o lbl.arts/ascii.o $DLFLAGS

./time_loop
