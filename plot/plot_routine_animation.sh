CC=gcc

CFLAGS='-O3 -ffast-math -fstrict-aliasing -g -Wall -fsanitize=address'

DLFLAGS='-lm'

$CC -c $CFLAGS -I./gnuplot/gnuplot_i-2.10/src/ -I./lbl.arts plot_routine_animation.c

$CC $CFLAGS -o plot_routine_animation plot_routine_animation.o gnuplot/gnuplot_i-2.10/gnuplot_i.o lbl.arts/ascii.o $DLFLAGS

./plot_routine_animation
