#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "vecmat.h"
#include "levy_fem.h"
#include "chi.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>

#include <omp.h>

#include "europeanopt.h"
#include "americanopt.h"
#include "swing.h"
#include "timer.h"

using namespace flens;
using namespace std;

// time-dependent jump intensity
double rho(double t){
	//return (2.0/(abs(sin(M_PI*(365*(1-t) -91)/365)) +1) -1);
	return 1;
} 

int
main(int argc, char *argv[])
{

	if (argc != 8){
		cout << "usage: " << argv[0] << " n m p R G M Y" << endl << endl
		     << "n: discretization level of the log price domain" << endl
		     << "m: discretization level of the time interval" << endl
		     << "R: truncation parameter" << endl
		     << "G,M,Y: parameters of the (CGMY) jump intensity" << endl;
			
		exit(EXIT_FAILURE);
	}


	int	n=atoi(argv[1]),	// discretization level of the log spot price
		m=atoi(argv[2]),	// discretization level of the time interval
		p=atoi(argv[3]);	// number of exercise rights
	
	

	double	sigma = 0.5,		// volatility
		lambda=0.2,		// mean reversion speed
		r = 0.03,		// interest rate
		tau = 1,		// time to maturity
		R = atof(argv[4]),	// truncation parameter
		a=-R, b=R,
		theta = 0.5,
		C_cgmy=1,
		G_cgmy=atof(argv[5]),
		M_cgmy=atof(argv[6]),
		Y_cgmy=atof(argv[7]),
		dt=tau/m,
		delta=((int) (0.2*tau/dt))*dt, // refraction period
		h=(exp(b)-exp(a))/(n-1),
		//h=(b-a)/(n-1),
		
		gamma=C_cgmy*gsl_sf_gamma(-Y_cgmy)* ( std::pow(M_cgmy-1,Y_cgmy)
		+std::pow(1+G_cgmy,Y_cgmy) -std::pow(G_cgmy,Y_cgmy) -std::pow(M_cgmy,Y_cgmy) )
		+(-C_cgmy*chi(1,Y_cgmy,G_cgmy) + C_cgmy*chi(1,Y_cgmy,M_cgmy)) - (0.5*sigma*sigma-r);
	
	char filename[100]="";
		

	DEVector::IndexVariable i;
	DEVector x(n),  tgrid(m+1);
	
	//set number of threads
	omp_set_num_threads(OMPTHREADS);
	
	// construct time grid
	tgrid(i) = (i-1)*dt;
	
	// construct space grid
	//x(i) = a + (i-1)*h;
	for (int j=1; j <= x.length(); j++)
		x(j) = log(exp(a) + (j-1)*h);



	cout << "gamma: " << gamma << endl;
	cout  << n << "  " << m << "  "  << R << "  " \
	<< G_cgmy << "  " << M_cgmy << "  " << Y_cgmy << endl;
	
	
	/* Configuration of the payoff function:
	 * The overall payoff function consits of the linear combination
	 * of put and call payoffs with strikes defined in put_strikes
	 * and call_strikes. Negative signs of strikes result in 
	 * weights of -1 of the corresponding payoff function.
	 * Example: 
	 * put_strikes={10, -20};
	 * call_strikes={50};
	 * Resulting payoff: max(10-S,0) - max(20-S,0) + max(S-50,0)
	 */
	vector<double> put_strikes = {100};//101, -50.5
	vector<double> call_strikes ={};//151.5, -202
	

	
	Swing option(x,tgrid,p,delta,r,sigma,gamma,lambda,theta,rho,G_cgmy,M_cgmy, Y_cgmy,put_strikes,call_strikes);
	//EuropeanOpt option(x,tgrid,r,sigma,gamma,lambda,theta,rho,G_cgmy,M_cgmy,Y_cgmy,put_strikes,call_strikes);
	//AmericanOpt option(x,tgrid,r,sigma,gamma,lambda,rho,G_cgmy,M_cgmy,Y_cgmy,put_strikes,call_strikes);
	

	option.initialize();
	option.solve();

	snprintf(filename, 100, "val_n%d_m%d_R%f_G%f_M%f_Y%f", n,m,R,G_cgmy,M_cgmy,Y_cgmy);
	option.write_to_file(filename);
	//option.write_to_file("val");

		
	return 0;
	
}



