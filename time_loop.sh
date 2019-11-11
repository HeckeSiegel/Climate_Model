CC=gcc

CFLAGS='-O3 -ffast-math -fstrict-aliasing -g -Wall'

DLFLAGS='-lm'

$CC -c $CFLAGS thermal_radiation.c


$CC -c $CFLAGS -I./gnuplot/gnuplot_i-2.10/src/ time_loop.c

$CC $CFLAGS -o time_loop time_loop.o thermal_radiation.o gnuplot/gnuplot_i-2.10/gnuplot_i.o $DLFLAGS

./time_loop
