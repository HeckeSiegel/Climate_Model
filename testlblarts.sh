CC=gcc

CFLAGS='-O3 -ffast-math -fstrict-aliasing -g -Wall'

DLFLAGS='-lm'

$CC -c $CFLAGS -I./lbl.arts testlblarts.c

$CC $CFLAGS -o testlblarts testlblarts.o lbl.arts/ascii.o $DLFLAGS

./testlblarts
