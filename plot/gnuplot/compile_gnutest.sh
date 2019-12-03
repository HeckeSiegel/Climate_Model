gcc -c -Wall  -I./gnuplot_i-2.10/src/ -g gnutest.c
gcc -o gnutest gnutest.o gnuplot_i-2.10/gnuplot_i.o 
