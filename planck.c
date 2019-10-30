#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define h 6.626e-34 //Js
#define c 2.99e8 // m/s
#define kB 1.38e-23 // m^2kgs^-2K^-1
#define sigma 5.67e-8
#define pi 3.14

double B_gray(double T){
    return sigma*pow(T,4)/pi;
}
double B_int(double lambda_1, double lambda_2, double T){
    int N = 100;
    double B = 0;
    double dwvl = (lambda_2 - lambda_1)/(N+1);
    double wvl[N];
    //initialize wavelengths for integration
    for(int i = 0; i<N; i++){
        wvl[i] = lambda_1 + (i+1)*dwvl;
    }
    //integration
    for(int i=0; i<N; i++){
        B += (2*h*pow(c,2)/pow(wvl[i],5) * 1./(exp(h*c/(wvl[i]*kB*T))-1.))*dwvl;
    }
    return B;
}

