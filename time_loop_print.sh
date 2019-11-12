CC=gcc

CFLAGS='-O3 -ffast-math -fstrict-aliasing -g -Wall'

DLFLAGS='-lm'

$CC -c $CFLAGS thermal_radiation.c


$CC -c $CFLAGS time_loop_print.c

$CC $CFLAGS -o time_loop_print time_loop_print.o thermal_radiation.o $DLFLAGS

./time_loop_print
