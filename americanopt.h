#ifndef AMERICANOPT_H
#define AMERICANOPT_H
#include "vecmat.h"
#include "levy_fem.h"
#include "derivative.h"
#include "psor.h"
#include <vector>
#include <string>
#include <iostream>
#include "signum.h"
#include <cmath>

class AmericanOpt: public Derivative{

private:
	size_t n, m;
	GEMatrix D, M,  Cl, J, Z;
	double K, r, sig, gamma, lambda, 
		 /*C_cgmy,*/ G_cgmy, M_cgmy, Y_cgmy;
	Func1D C_cgmy;
	vector<double> put_strikes, call_strikes;


public:


	AmericanOpt(DEVector grid_, DEVector tgrid_, 
	double r_, double sig_, double gamma_, double lambda_, 
	Func1D C_, double G_, double M_, double Y_,
	vector<double> put_strikes_, vector<double> call_strikes_)
	: Derivative(grid_, tgrid_,std::vector<DEVector>(tgrid_.length(),DEVector(grid_.length()))),
	 n(grid_.length()), m(tgrid_.length()-1),
	 D(GEMatrix(n,n)), M(GEMatrix(n,n)), J(GEMatrix(n,n)),
	 r(r_), sig(sig_), gamma(gamma_), lambda(lambda_),
	 G_cgmy(G_), M_cgmy(M_), Y_cgmy(Y_), C_cgmy(C_),
	 put_strikes(put_strikes_), call_strikes(call_strikes_) {}
	 
	 void initialize();
	 void solve();
	 void write_to_file(std::string filename);
	
};

void AmericanOpt::initialize(){
	// assemble mass and stiffness matrices
	assemble_mass_mtrx(grid,  M);
	assemble_diff_op(sig, gamma, lambda,  grid,  D);
	//Assemble jump matrix with C=1. => Z=D + C(t)*J-r*M;
	assemble_jump_mat(grid, J, 1,G_cgmy, M_cgmy, Y_cgmy);
	//Z=D + J-r*M;
}

// solve the sequence of P(I)DIs to get the price of the American option
void AmericanOpt::solve(){

	double	dt, theta=1, rho;
	DEVector rhs(n), F(n), c0(n),tmp(n);
	GEMatrix Cl(n,n),Cr(n,n);
	int it;

	c0=0;
	

	
	for (size_t j=1; j <= m; j++){
		dt= tgrid(j+1) -tgrid(j);
		rho = C_cgmy(j*dt); // rho = timedependent C_cgmy
		Z=D + rho*J-r*M;
	
		
		Cr=(M+dt*(1.0-theta)*Z);
		Cl=(M-dt*theta*Z);
		
		
		// assemble rhs
		F=0;
		for(vector<double>::iterator it = put_strikes.begin(); it != put_strikes.end(); ++it){
			assemble_rhs_am_put(rhs, grid, 0.0, sig, gamma, 
				lambda, r,  fabs(*it),  rho,  G_cgmy,  M_cgmy,  Y_cgmy);
  			F += signum(*it)*rhs;
		}
		for(vector<double>::iterator it = call_strikes.begin(); it != call_strikes.end(); ++it){
			assemble_rhs_am_call(rhs, grid, 0.0, sig, gamma,
				 lambda, r,  fabs(*it),  rho,  G_cgmy,  M_cgmy,  Y_cgmy);
  			F += signum(*it)*rhs;
		}


		rhs = Cr*u[j-1] +dt*F ;


		it=psor(Cl,u[j],u[j-1],rhs,c0);

		cout << "iteration: " << j << " / " << m << " iterations: " << it << endl;
	}
}

// calculate the payoff at point x
double am_payoff(double x, vector<double> put_strikes, vector<double> call_strikes){
	double ret=0;


	for(vector<double>::iterator it = put_strikes.begin(); it != put_strikes.end(); ++it){
  		ret += signum(*it)*std::max(fabs(*it)-exp(x),0.0);
	}
	for(vector<double>::iterator it = call_strikes.begin(); it != call_strikes.end(); ++it){
  		ret += signum(*it)*std::max(exp(x) -fabs(*it),0.0);
	}
	return ret;

}

// write results to output file
void AmericanOpt::write_to_file(std::string filename){
	ofstream file(filename);
	for (int j=1; j <= u[m].length(); j++){
		file << exp(grid(j)) << "  " << u[m](j) + am_payoff(grid(j),put_strikes, call_strikes)  
		<< "  " << am_payoff(grid(j),put_strikes, call_strikes) << endl;
	}
}

#endif
