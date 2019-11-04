CC=gcc

CFLAGS='-O2 -g -Wall'

DLFLAGS='-lm'

$CC -c $CFLAGS thermal_radiation.c


$CC -c $CFLAGS thermal_radiation_test.c

$CC $CFLAGS -o thermal_radiation_test thermal_radiation_test.o thermal_radiation.o $DLFLAGS

./thermal_radiation_test
