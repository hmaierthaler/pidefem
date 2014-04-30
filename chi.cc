#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <iostream>

double chi(double x, double alpha, double L){

	return (pow(L,alpha-1.)*gsl_sf_gamma_inc(3.0-alpha, L*x)
		+(-2.+alpha)*exp(-L*x)*pow(x,-alpha+1.)
		-L*exp(-L*x)*pow(x,2.-alpha))/(2.+alpha*alpha-3.*alpha);
}


double hat_func(double x, double x0, double x1, double x2){

	if (x < x0 || x > x2){
		return 0.0;
	}

	if (x < x1){
		return (x-x0)/(x1-x0);
	}else{
		return (x2-x)/(x2-x1);
	}
}
