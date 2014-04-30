#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf.h>
#include "chi.h"
#include "vecmat.h"
#include <omp.h>

typedef double (*Func1D)(double);

// needed for Levy part
double I_rhs_call(double x, void *param){
	double *alpha	= (double*) param,
	       tau	= alpha[0],
	       //sigma	= alpha[1],
	       r	= alpha[2],
	       K	= alpha[3],
	       C	= alpha[4],
	       G	= alpha[5],
	       M	= alpha[6],
	       Y	= alpha[7],

	       bu	= -x -r*tau + log(K),
	       ret	=0;

	if (bu < 0){

		ret =	-(exp(x)-exp(-r*tau)*K)*C*chi(-bu,1+Y,G) 
			+ exp(x)*C*chi(-bu,Y,G);
	} else {
		ret = 	exp(x)*C*chi(bu,1.0+Y,M-1.0) 
			-exp(-r*tau)*K*C*chi(bu,1.0+Y,M);
		      
	}

	ret *= alpha[8] + x*alpha[9];
	return ret;

}

// assemble the load vector for a European call option
void assemble_rhs_eur_call(	DEVector &rhs, DEVector grid, double offset,
		 		double tau, double sigma, double gamma, 
		 		double lambda, double r, double K,
		 		double C, double G, double M, double Y){
		 	
	gsl_integration_workspace *w[OMPTHREADS];
	gsl_function F;
	double x0,x1,x2,
	       res[OMPTHREADS],
	       err[OMPTHREADS],
	       *param[OMPTHREADS],
	       *sing[OMPTHREADS],
	       bu = -r*tau+log(K) -offset,
	       h;

	int n=rhs.length(),
	    id;
		
	size_t	si = 4000;


	for (int i=0; i < OMPTHREADS; i++){
		w[i] = gsl_integration_workspace_alloc(4000);
		param[i] = new double[10];
		sing[i] = new double[2];
	}
	

	#pragma omp parallel for private(id,F,x0,x1,x2,h)
	for (int i=1; i<=n; i++){
	
		id = omp_get_thread_num();
		
		param[id][0] = tau;
		param[id][1] = sigma;
		param[id][2] = r;
		param[id][3] = K;
		param[id][4] = C;
		param[id][5] = G;
		param[id][6] = M;
		param[id][7] = Y;
	
	
		// Differential part
		x0=(i>1) ? grid(i-1) : grid(i);
		x1=grid(i);
		x2=(i<n) ? grid(i+1) : grid(i);
		rhs(i)=+0.5*sigma*sigma*K*exp(-r*tau)*hat_func(bu ,x0,x1,x2);
		
		if (i>1 && i<n && bu < x0){
			rhs(i) += (gamma+0.5*sigma*sigma-r)*((x2-x1)*exp(x0)
				+(x0-x2)*exp(x1)-exp(x2)*(x0-x1))/((x0-x1)*(x1-x2))
				-lambda*(-(x1-x2)*(x0-2)*exp(x0)
				+(x1-2)*(x0-x2)*exp(x1)-exp(x2)*(x2-2)*(x0-x1))/((x0-x1)*(x1-x2));
			
		}
		
		if (i > 1 && x0 < bu && bu < x1){
			rhs(i) += (gamma+0.5*sigma*sigma-r)*((x1-x2)*(bu-x0-1)*exp(bu)
				+(x0-x2)*exp(x1)-exp(x2)*(x0-x1))/((x0-x1)*(x1-x2)) 
				-lambda*((x1-x2)*((-bu+1)*x0+bu*bu-2*bu+2)*exp(bu)
				+(x1-2)*(x0-x2)*exp(x1)
				-exp(x2)*(x2-2)*(x0-x1))/((x0-x1)*(x1-x2));
		}

		if (i < n && x1 < bu && bu < x2){
			rhs(i) += (gamma+0.5*sigma*sigma-r)*((-bu+x2+1)*exp(bu)-exp(x2))/(x1-x2)
				-lambda*(((bu-1)*x2-bu*bu+2*bu-2)*exp(bu)-exp(x2)*(x2-2))/(x1-x2);
		}


		

		//// Levy part
		// integral from x_{i‐1} to x_{i}
		F.function = &I_rhs_call;
		if (i>1){
			param[id][8] = -grid(i-1)/(grid(i)-grid(i-1));
			param[id][9] = 1.0/(grid(i)-grid(i-1));
			F.params=param[id];
			h=grid(i)-grid(i-1);

				sing[id][0] = grid(i-1);
				sing[id][1] = grid(i);

				gsl_integration_qags(&F,sing[id][0]+1e-7*h,
					sing[id][1]-1e-7*h,0,1e-4,si, w[id],&res[id],&err[id]);
			rhs(i) += res[id];

		}
			
		// integral from x_{i‐1} to x_{i}
		if (i<n){
			param[id][8] = grid(i+1)/(grid(i+1)-grid(i));
			param[id][9] = -1.0/(grid(i+1)-grid(i));
			F.params=param[id];
			
			h=grid(i+1)-grid(i);

				sing[id][0] = grid(i);
				sing[id][1] = grid(i+1);

				gsl_integration_qags(&F,sing[id][0]+1e-7*h,
					sing[id][1]-1e-7*h,0,1e-4,si, w[id],&res[id],&err[id]);

			rhs(i) +=res[id];

		}

	}
	
	for (int i=0; i < OMPTHREADS; i++){
		gsl_integration_workspace_free(w[i]);
		delete [] sing[i];
		delete [] param[i];
	}
}


// assemble the load vector for an American call option
void assemble_rhs_am_call(DEVector &rhs, DEVector grid, double offset,
			double sigma, double gamma, double lambda, double r, double K,
			double C, double G, double M, double Y){



	gsl_integration_workspace *w[OMPTHREADS];
	gsl_function F;
	double x0,x1,x2,
	       res[OMPTHREADS],
	       err[OMPTHREADS],
	       *param[OMPTHREADS],
	       *sing[OMPTHREADS],
	       bu = log(K) -offset,
	       h;

	int n=rhs.length(),
	    id;
	size_t	si = 4000;

	for (int i=0; i < OMPTHREADS; i++){
		w[i] = gsl_integration_workspace_alloc(4000);
		param[i] = new double[10];
		sing[i] = new double[2];
	}
	

	
	#pragma omp parallel for private(id,F,x0,x1,x2,h)
	for (int i=1; i<=n; i++){

		id = omp_get_thread_num();
		
		// Differential part
		x0=(i>1) ? grid(i-1) : grid(i);
		x1=grid(i);
		x2=(i<n) ? grid(i+1) : grid(i);
		rhs(i)=+0.5*sigma*sigma*K*hat_func(bu ,x0,x1,x2);
		
		if (i>1 && i<n && bu < x0){
			rhs(i) += (gamma+0.5*sigma*sigma-r)*((x2-x1)*exp(x0)
				+(x0-x2)*exp(x1)-exp(x2)*(x0-x1))/((x0-x1)*(x1-x2))
				-lambda*(-(x1-x2)*(x0-2)*exp(x0)+(x1-2)*(x0-x2)*exp(x1)
				-exp(x2)*(x2-2)*(x0-x1))/((x0-x1)*(x1-x2))
				-r*K*(-0.5*x0+0.5*x2);
			
		}
		
		if (x0 < bu && bu < x1){
			rhs(i) += (gamma+0.5*sigma*sigma-r)*((x1-x2)*(bu-x0-1)*exp(bu)
				+(x0-x2)*exp(x1)-exp(x2)*(x0-x1))/((x0-x1)*(x1-x2)) 
				-lambda*((x1-x2)*((-bu+1)*x0+bu*bu-2*bu+2)*exp(bu)
				+(x1-2)*(x0-x2)*exp(x1)
				-exp(x2)*(x2-2)*(x0-x1))/((x0-x1)*(x1-x2))
				+r*K*((-2*bu+x1+x2)*x0+bu*bu-x2*x1)/(2*x0-2*x1);
		}

		if (x1 < bu && bu < x2){
			rhs(i) += (gamma+0.5*sigma*sigma-r)*((-bu+x2+1)*exp(bu)-exp(x2))/(x1-x2) 
				-lambda*(((bu-1)*x2-bu*bu+2*bu-2)*exp(bu)-exp(x2)*(x2-2))/(x1-x2)
				-r*K*(bu-x2)*(bu-x2)/(2*x1-2*x2);
		}
		


		//// Jump part
		param[id][0] = 0.0; //tau
		param[id][1] = sigma;
		param[id][2] = r; // r 
		param[id][3] = K;
		param[id][4] = C;
		param[id][5] = G;
		param[id][6] = M;
		param[id][7] = Y;
		F.function = &I_rhs_call;
		
		// integral from x_{i‐1} to x_{i}
		if (i>1){
			param[id][8] = -grid(i-1)/(grid(i)-grid(i-1));
			param[id][9] = 1.0/(grid(i)-grid(i-1));
			F.params=param[id];
			h=grid(i)-grid(i-1);

				sing[id][0] = grid(i-1);
				sing[id][1] = grid(i);

				gsl_integration_qags(&F,sing[id][0]+1e-6*h,
					sing[id][1]-1e-6*h,0,1e-4,si, w[id],&res[id],&err[id]);

			rhs(i) += res[id];

		}
			
		// integral from x_{i} to x_{i+1}
		if (i<n){
			param[id][8] = grid(i+1)/(grid(i+1)-grid(i));
			param[id][9] = -1.0/(grid(i+1)-grid(i));
			F.params=param[id];
			
			h=grid(i+1)-grid(i);
	
				sing[id][0] = grid(i);
				sing[id][1] = grid(i+1);

				gsl_integration_qags(&F,sing[id][0]+1e-6*h,
					sing[id][1]-1e-6*h,0,1e-4,si, w[id],&res[id],&err[id]);

			rhs(i) +=res[id];


		}

	}
	//gsl_integration_workspace_free(w);
	for (int i=0; i < OMPTHREADS; i++){
		gsl_integration_workspace_free(w[i]);
		delete [] sing[i];
		delete [] param[i];
	}

}

// assemble the load vector for the early exercise part of the swing option
void assemble_rhs_swing_am_call(DEVector &rhs, DEVector grid, double offset,
				double tau, double delta, int p, double sigma, 
				double gamma, double lambda, double r, double K,
				Func1D C, double G, double M, double Y){
				
	DEVector rhs2(rhs.length());
	
	if (p==0){
		rhs=0;
		return;
	}
		
	
	// p exercises possible
	if ( tau >= (p-1.0)*delta){

			assemble_rhs_am_call(rhs, grid, 0.0, sigma, 
				gamma, lambda, r, K, C(tau), G, M, Y);

			assemble_rhs_swing_am_call(rhs2, grid, r*delta, tau-delta, 
				delta, p-1,sigma,gamma,lambda,r,K,C,G,M, Y);
			rhs = rhs + exp(-r*delta) *rhs2;
	}
	// p-1 exercises possible
	else if (tau >=0 && tau < (p-1.0)*delta) {
		assemble_rhs_swing_am_call(rhs, grid, 0.0, tau, delta, p-1, 
				sigma, gamma, lambda, r, K, C,G,M,Y);
	}
	

}

// assemble the load vector for the "European part" of the swing option
void assemble_rhs_swing_eur_call(DEVector &rhs, DEVector grid, double offset,
				double tau, double delta, int p, double sigma, 
				double gamma, double lambda, double r, double K,
				Func1D C, double G, double M, double Y){
				
	DEVector rhs1(rhs.length()), rhs2(rhs.length());
	
	if (p==0){
		rhs=0;
		return;
	}
		
	
	// p exercises possible
	if ( tau >= (p-1-0)*delta){
			// 6th parameter = gamma
			assemble_rhs_eur_call(rhs1, grid, 0.0, tau, 
				sigma, gamma, lambda, r, K, C(tau), G, M, Y);

			assemble_rhs_swing_eur_call(rhs2, grid, r*delta, tau,
				delta, p-1, sigma, gamma, lambda, r, K, C, G, M, Y);
			rhs = rhs1 + exp(-r*delta) *rhs2;
	}
	// p-1 exercises possible
	else if (tau >=0 && tau < (p-1.0)*delta) {
		assemble_rhs_swing_eur_call(rhs, grid, 0.0, tau, delta, p-1, 
			sigma, gamma, lambda, r, K, C,G,M,Y);
	}
	

}


