#ifndef SWING_H
#define SWING_H
#include "vecmat.h"
#include "levy_fem.h"
#include "derivative.h"
#include "psor.h"
#include "signum.h"
#include <vector>
#include <string>
#include <iostream>

class Swing: public Derivative{

private:
	size_t n, m;
	int p;
	GEMatrix D, M,   J, Z;
	double delta, r, sig, gamma, lambda, theta,
		 /*C_cgmy,*/ G_cgmy, M_cgmy, Y_cgmy;
	Func1D C_cgmy;
	vector<double> put_strikes, call_strikes;


public:


	Swing(DEVector grid_, DEVector tgrid_, 
	int p_, double delta_,
	double r_, double sig_, double gamma_, double lambda_, double theta_,
	Func1D C_, double G_, double M_, double Y_, vector<double> put_strikes_, vector<double> call_strikes_)
	: Derivative(grid_, tgrid_,std::vector<DEVector>(tgrid_.length(),DEVector(grid_.length()))),
	 n(grid_.length()), m(tgrid_.length()-1), p(p_), 
	 D(GEMatrix(n,n)), M(GEMatrix(n,n)), J(GEMatrix(n,n)),
	 delta(delta_), r(r_), sig(sig_), gamma(gamma_), lambda(lambda_), theta(theta_),
	 G_cgmy(G_), M_cgmy(M_), Y_cgmy(Y_),C_cgmy(C_),
	 put_strikes(put_strikes_), call_strikes(call_strikes_) {}
	 
	void initialize();
	void solve();
	void write_to_file(std::string filename);
	void write_ex_strat(std::string filename);
};

void Swing::initialize(){
	// assemble mass and stiffness matrices
	assemble_mass_mtrx(grid,  M);
	assemble_diff_op(sig, gamma, lambda,  grid,  D);
	assemble_jump_mat(grid, J, 1,G_cgmy, M_cgmy, Y_cgmy);

	//J=0;
	//Z=D + J-r*M;
}

void Swing::solve(){

	double	dt=tgrid(2) -tgrid(1),rho;
	DEVector rhs(n), F(n), c0(n),tmp(n),rhs2(n);
	GEMatrix Cr, Cl ,Am;
	IndexVector piv(n);
	

	int 	it=0,
		ndelta = delta/dt;

	// implementation of Algorithm 2
	for (int l=1; l <= p; l++){
		for ( size_t j = (l-1)*ndelta; j <= m; j++){
			if (j==0){
				continue;
			}
			dt=tgrid(j+1)-tgrid(j);
			rho=C_cgmy(dt*j);
			Z=D + rho*J-r*M;
			Am = M-dt*Z;
			Cr=(M+dt*(1.0-theta)*Z);
			Am = M-dt*Z;
				
			if (l > 1) {
				// European option pricing problem
				// solve (5.66)
				
				// assemble load vector
				rhs=0;
				for(vector<double>::iterator it = put_strikes.begin(); it != put_strikes.end(); ++it){
  					assemble_rhs_swing_eur_put(rhs2,grid, 0,(j+1)*dt ,delta,l-1,
  							sig,gamma,lambda,r,fabs(*it),C_cgmy,G_cgmy,M_cgmy,Y_cgmy);
  					rhs += signum(*it)*rhs2;
				}
				for(vector<double>::iterator it = call_strikes.begin(); it != call_strikes.end(); ++it){
  					assemble_rhs_swing_eur_call(rhs2,grid, 0,(j+1)*dt,delta,l-1,
  							sig,gamma,lambda,r,fabs(*it),C_cgmy,G_cgmy,M_cgmy,Y_cgmy);
  					rhs += signum(*it)*rhs2;
				}

				c0 = Cr*u[j-ndelta] +dt*rhs;
				
				for (int k=1; k <= ndelta; k++){
					if (k>1){
						// assemble load vector
						rhs=0;
						for(vector<double>::iterator it = put_strikes.begin(); it != put_strikes.end(); ++it){
							rhs2=0;
  							assemble_rhs_swing_eur_put(rhs2, grid, 0,(j+k)*dt,delta,
									l-1,sig,gamma,lambda, r,fabs(*it),C_cgmy,G_cgmy,M_cgmy,Y_cgmy);
  							rhs += signum(*it)*rhs2;
						}
						for(vector<double>::iterator it = call_strikes.begin(); it != call_strikes.end(); ++it){
							rhs2=0;
  							assemble_rhs_swing_eur_call(rhs2, grid, 0,(j+k)*dt,delta,
									l-1,sig,gamma,lambda, r,fabs(*it),C_cgmy,G_cgmy,M_cgmy,Y_cgmy);
  							rhs += signum(*it)*rhs2;
						}
						tmp = Cr*c0;
						c0 = tmp +dt*rhs;
					}
					// lapack::sv manipulates the Matrix!!!
					Cl=(M-dt*theta*Z);
					lapack::sv(Cl,piv,c0);
				}
				
			} else {
				c0 = 0;
			}
			
			if (j>0) {
				// American option pricing problem
				// solve 5.65
				rhs=0;
				for(vector<double>::iterator it = put_strikes.begin(); it != put_strikes.end(); ++it){
					rhs2=0;
  					assemble_rhs_swing_am_put(rhs2,grid,0,(j)*dt,delta,l,
  							sig,gamma,lambda,r,fabs(*it),C_cgmy,G_cgmy,M_cgmy,Y_cgmy);
  					rhs+=signum(*it)*rhs2;
				}
				for(vector<double>::iterator it = call_strikes.begin(); it != call_strikes.end(); ++it){
					rhs2=0;
  					assemble_rhs_swing_am_call(rhs2,grid,0,(j)*dt,delta,l,
  							sig,gamma,lambda,r,fabs(*it),C_cgmy,G_cgmy,M_cgmy,Y_cgmy);
  					rhs+=signum(*it)*rhs2;
				}
				
				Cl=(M-dt*theta*Z);
				F = Cr*u[j-1] +dt*rhs ;
				it=psor(Cl,u[j],c0,F,c0);
				cout << " l=" << l << "  j=" << j<< " it: " << it << endl;
			}
		}	 
	}
}


double multi_payoff(	double x, double tau, int p, double delta, double r, 
			vector<double> put_strikes, vector<double> call_strikes){
	double ret=0;
	
	if (p==0){
		return 0.0;
	}
	
	if (tau >= (p-1.0)*delta){
		for(vector<double>::iterator it = put_strikes.begin(); it != put_strikes.end(); ++it){
  			ret += signum(*it)*std::max(fabs(*it)-exp(x),0.0);
		}
		for(vector<double>::iterator it = call_strikes.begin(); it != call_strikes.end(); ++it){
  			ret += signum(*it)*std::max(exp(x) -fabs(*it),0.0);
		}
			
		return ( ret + exp(-r*delta)*multi_payoff(x+r*delta, tau-delta, p-1, 
			delta, r, put_strikes, call_strikes) );
	} else {
		return multi_payoff( x,  tau,  p-1,  delta,  r, put_strikes, call_strikes);
	}
}

void Swing::write_to_file(std::string filename){
	ofstream file(filename);
	
	
	for (size_t j=1; j<=m+1;j++){
		for (size_t k=1; k <=n; k++){
			file << tgrid(j) << "  "<< exp(grid(k)) << "  " \
			<< u[j-1](k) + multi_payoff(grid(k), tgrid(j), p, delta, r, put_strikes, call_strikes) << endl;
		}
		file << endl;
	}
	
	/*
	double uz=0;
	for (size_t k=1; k<=n; k++){
		
		uz=0;
		for(vector<double>::iterator it = put_strikes.begin(); it != put_strikes.end(); ++it){
  			uz += signum(*it)*std::max(fabs(*it)-exp(grid(k)),0.0);
		}
		for(vector<double>::iterator it = call_strikes.begin(); it != call_strikes.end(); ++it){
  			uz += signum(*it)*std::max(exp(grid(k)) -fabs(*it),0.0);
		}
		
		
		file << exp(grid(k)) << "  " << u[m](k) + multi_payoff(grid(k), tgrid(m+1), p, delta, r, put_strikes, call_strikes)
		<< "  " << uz << endl;
	}*/
	
	/*
	for (size_t j=1; j <=m+1;j++){
		file << tgrid(j)<< "  "<< u[j-1](69) << endl;
	}*/
}

#endif
