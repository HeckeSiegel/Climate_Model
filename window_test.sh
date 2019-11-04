CC=gcc

CFLAGS='-O3 -ffast-math -g -Wall'

DLFLAGS='-lm'

$CC -c $CFLAGS thermal_radiation.c


$CC -c $CFLAGS window_test.c

$CC $CFLAGS -o window_test window_test.o thermal_radiation.o $DLFLAGS

./window_test
