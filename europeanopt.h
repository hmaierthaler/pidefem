#ifndef EUROPEANOPT_H
#define EUROPEANOPT_H
#include "vecmat.h"
#include "levy_fem.h"
#include "derivative.h"
#include <vector>
#include <string>
#include <iostream>
#include "signum.h"
#include "sor.h"
#include <cmath>

class EuropeanOpt: public Derivative{

private:
	size_t n, m;
	GEMatrix D, M, J, Z;
	double r, sig, gamma, lambda, theta,
		 /*C_cgmy, */G_cgmy, M_cgmy, Y_cgmy;
	Func1D C_cgmy;
	
	vector<double> put_strikes,  call_strikes;


public:


	EuropeanOpt(DEVector grid_, DEVector tgrid_, 
	double r_, double sig_, double gamma_, double lambda_, double theta_,
	/*double C_,*/ Func1D C_, double G_, double M_, double Y_,
	vector<double> put_strikes_, vector<double> call_strikes_)
	: Derivative(grid_, tgrid_,std::vector<DEVector>(tgrid_.length(),DEVector(grid_.length()))),
	 n(grid_.length()), m(tgrid_.length()-1),
	 D(GEMatrix(n,n)), M(GEMatrix(n,n)), J(GEMatrix(n,n)),
	 r(r_), sig(sig_), gamma(gamma_), lambda(lambda_), theta(theta_),
	 G_cgmy(G_), M_cgmy(M_), Y_cgmy(Y_),C_cgmy(C_),
	 put_strikes(put_strikes_), call_strikes(call_strikes_) {}
	 
	 void initialize();
	 void solve();
	 void write_to_file(std::string filename);
	
};

void EuropeanOpt::initialize(){
	// assemble mass and stiffness matrix
	assemble_mass_mtrx(grid,  M);
	assemble_diff_op(sig, gamma, lambda,  grid,  D);
	//Assemble jump matrix with C=1. => Z=D + C(t)*J-r*M;
	assemble_jump_mat(grid, J, 1 ,G_cgmy, M_cgmy, Y_cgmy); 
	//Z=D + J-r*M;
}

void EuropeanOpt::solve(){

	double	dt,rho;
	DEVector F(n),rhs(n);
	GEMatrix Cr, Cl;
	IndexVector piv(n);
	//u=std::vector<DEVector>(m+1,DEVector(n));


	for (size_t j=1; j <= m; j++){
		cout << "iteration: " << j << " / " << m << endl;
		dt= tgrid(j+1) -tgrid(j);
		rho = C_cgmy(j*dt);
		Z=D + rho*J-r*M;
		
		Cr=(M+dt*(1.0-theta)*Z);
		Cl=(M-dt*theta*Z);

		
		F=0;
		for(vector<double>::iterator it = put_strikes.begin(); it != put_strikes.end(); ++it){
			assemble_rhs_eur_put(rhs,  grid,  0, tgrid(j),  sig, gamma ,  lambda,  r,  fabs(*it),
		 			rho,  G_cgmy,  M_cgmy,  Y_cgmy);
  			F += signum(*it)*rhs;
		}
		for(vector<double>::iterator it = call_strikes.begin(); it != call_strikes.end(); ++it){
			assemble_rhs_eur_call(rhs,  grid,  0, tgrid(j),  sig, gamma ,  lambda,  r,  fabs(*it),
		 			rho,  G_cgmy,  M_cgmy,  Y_cgmy);
  			F += signum(*it)*rhs;
		}

		u[j]=Cr*u[j-1] +dt*F ;

		// lapack::sv manipulates the Matrix!!!
		lapack::sv(Cl,piv,u[j]);
		//sor(Cl,u[j],u[j-1],rhs);
	}

}

double eur_payoff(double x, vector<double> put_strikes, vector<double> call_strikes){
	double ret=0;


	for(vector<double>::iterator it = put_strikes.begin(); it != put_strikes.end(); ++it){
  		ret += signum(*it)*std::max(fabs(*it)-exp(x),0.0);
	}
	for(vector<double>::iterator it = call_strikes.begin(); it != call_strikes.end(); ++it){
  		ret += signum(*it)*std::max(exp(x) -fabs(*it),0.0);
	}
	return ret;

}

void EuropeanOpt::write_to_file(std::string filename){
	ofstream file(filename);
	for (int j=1; j <= u[m].length(); j++){
		file << exp(grid(j)) << "  " << u[m](j) + exp(-r*tgrid(m))*eur_payoff(grid(j)+r*tgrid(m),put_strikes,call_strikes) 
		<< "  " <<  eur_payoff(grid(j),put_strikes,call_strikes) <<endl;
	}
}


#endif
