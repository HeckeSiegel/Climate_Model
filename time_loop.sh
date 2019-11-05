CC=gcc

CFLAGS='-O3 -ffast-math -fstrict-aliasing -g -Wall'

DLFLAGS='-lm'

$CC -c $CFLAGS thermal_radiation.c


$CC -c $CFLAGS time_loop.c

$CC $CFLAGS -o time_loop time_loop.o thermal_radiation.o $DLFLAGS

./time_loop
