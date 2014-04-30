#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf.h>
#include "vecmat.h"
#include "chi.h"
#include <omp.h>


double I_func(double x, void *params){

	double *tmp = (double *) params;

	double	x0 = tmp[0],
		x1 = tmp[1],
		x2 = tmp[2],
		h1 = x1-x0,
		h2 = x2-x1,
		C = tmp[3],
		G = tmp[4],
		M = tmp[5],
		Y = tmp[6],
		dx0 = x-x0,
		dx1 = x-x1,
		dx2 = x-x2,
		res=0;

	// case 1
	if (x > x2){
		res =	C/h1*(-dx0*(chi(dx0, 1.0+Y,G) - chi(dx1, 1.0+Y, G)) 
			+ chi(dx0,Y,G) - chi(dx1,Y,G))
			+C/h2*(dx2*( chi(dx1, 1.0+Y, G) - chi(dx2, 1.0+Y, G) ) 
			- (chi(dx1, Y, G) -chi(dx2, Y, G)));
	}
	// case 2
	if (x> x1 && x < x2){
		res =	+C/h1*(-dx0*( chi(dx0, 1.0+Y, G) - chi(dx1, 1.0+Y, G) )
			 + chi(dx0, Y, G) -chi(dx1, Y, G))
			+C*dx2/h2*chi(dx1,1+Y,G) - C/h2*chi(dx1,Y,G)
			+C*dx2/h2*chi(-dx2,1+Y,M) + C/h2*chi(-dx2,Y,M);
	}
	// case 3
	if (x > x0 && x < x1){
		res =	-dx0/h1*C*chi(dx0,1+Y,G) + C/h1*chi(dx0,Y,G)
			+C/h2*(-dx2*(chi(-dx1,1.0+Y,M)-chi(-dx2,1.0+Y,M) ) 
			- (chi(-dx1,Y,M) - chi(-dx2,Y,M) ) )
			-dx0/h1*C*chi(-dx1,1+Y,M) -C/h1*chi(-dx1,Y,M);
			
	}
	// case 4
	if (x < x0){
		 res =	C/h1*( dx0*(chi(-dx0, 1.0+Y,M) - chi(-dx1,1.0+Y,M) )
		  	+chi(-dx0,Y,M)-chi(-dx1,Y,M))
		  	+C/h2*(-dx2*(chi(-dx1, 1.0+Y,M)-chi(-dx2,1.0+Y,M) ) 
		  	-(chi(-dx1,Y,M) -chi(-dx2,Y,M) ) );
	}

	res *= tmp[7] + x*tmp[8];
	return res;
}

// I(x) for j=1 (lower boundary)
double I_func_lbound(double x, void *params){

	double *tmp = (double *) params;

	double	x1 = tmp[1], 
		x2 = tmp[2],

		h2 = x2-x1,
		C = tmp[3],
		G = tmp[4],
		M = tmp[5],
		Y = tmp[6],

		dx1 = x-x1,
		dx2 = x-x2,
		res=0;

	// case 1
	if (x > x2){
		res =	+C/h2*(dx2*( chi(dx1, 1.0+Y, G) - chi(dx2, 1.0+Y, G) ) 
			- (chi(dx1, Y, G) -chi(dx2, Y, G)));
	}
	// case 2
	if (x> x1 && x < x2){
		res =	+C*dx2/h2*chi(dx1,1+Y,G) - C/h2*chi(dx1,Y,G) 
			+C*dx2/h2*chi(-dx2,1+Y,M) + C/h2*chi(-dx2,Y,M);
	}


	res *= tmp[7] + x*tmp[8];
	return res;
}

// I(x) for j=n (upper boundary)
double I_func_ubound(double x, void *params){

	double *tmp = (double *) params;

	double	x0 = tmp[0],
		x1 = tmp[1],
		//x2 = tmp[2],
		h1 = x1-x0,
		//h2 = x2-x1,
		C = tmp[3],
		G = tmp[4],
		M = tmp[5],
		Y = tmp[6],
		dx0 = x-x0,
		dx1 = x-x1,
		res=0;


	// case 3
	if (x > x0 && x < x1){

		res =	-dx0/h1*C*chi(dx0,1+Y,G) + C/h1*chi(dx0,Y,G) 
			-dx0/h1*C*chi(-dx1,1+Y,M) -C/h1*chi(-dx1,Y,M);
			
	}
	// case 4
	if (x < x0){
		 res =	C/h1*( dx0*(chi(-dx0, 1.0+Y,M) - chi(-dx1,1.0+Y,M) ) 
		 	+ chi(-dx0,Y,M) -chi(-dx1,Y,M) );
	}

	res *= tmp[7] + x*tmp[8];
	return res;
}

// assemble the jump part of the stiffness matrix
void assemble_jump_mat(DEVector grid, GEMatrix &J, double C,double G, double M, double Y){
	gsl_integration_workspace *w[OMPTHREADS];
	gsl_function F;
	double	res[OMPTHREADS],
		err[OMPTHREADS],
		*param[OMPTHREADS],
		*sing[OMPTHREADS],
		h;
	size_t	si = 4000;
	int	n = grid.length(),
		id;


	for (int i=0; i < OMPTHREADS; i++){
		w[i] = gsl_integration_workspace_alloc(4000);
		param[i] = new double[9];
		sing[i] = new double[2];
	}

	for (int i=1; i <= n; i++) {
	
		#pragma omp parallel for private(id,F,h)
		for (int j=1; j <= n; j++){
			//cout << "i=" << i<< "  j=" << j << endl;
			id = omp_get_thread_num();
			
			F.function = &I_func;
			if (j==1){
				F.function = &I_func_lbound;
			}
			if (j==n){
				F.function = &I_func_ubound;
			}

			// set grid points x_j-1, x_j, x_j+1
			param[id][0] = (j > 1) ? grid(j-1) : 0;
			param[id][1] = grid(j);
			param[id][2] = (j < n) ? grid(j+1) : 0;
			param[id][3] = C;
			param[id][4] = G;
			param[id][5] = M;
			param[id][6] = Y;
			
			// integral from x_{i-1} to x_{i}
			if (i>1){
				h=grid(i)-grid(i-1);
				param[id][7] = -grid(i-1)/(grid(i)-grid(i-1));
				param[id][8] = 1.0/(grid(i)-grid(i-1));
				F.params=param[id];
				sing[id][0] = grid(i-1);
				sing[id][1] = grid(i);
				//gsl_integration_qagp(&F,sing,2,0,0.001,si, w,&res,&err);
				//if ( Y < 1)
				gsl_integration_qags(&F,sing[id][0]+1e-6*h,
					sing[id][1]-1e-6*h,0,1e-4,si, w[id],&res[id],&err[id]);
				J(i,j) += res[id];
				//cout << res;
			}
			
			// integral from x_{i} to x_{i+1}
			if (i<n){
				h=grid(i+1)-grid(i);
				param[id][7] = grid(i+1)/(grid(i+1)-grid(i));
				param[id][8] = -1.0/(grid(i+1)-grid(i));
				F.params=param[id];
				sing[id][0] = grid(i);
				sing[id][1] = grid(i+1); 
				//gsl_integration_qagp(&F,sing,2,0,0.001,si, w,&res,&err);
				gsl_integration_qags(&F,sing[id][0]+1e-6*h,
					sing[id][1]-1e-6*h,0,1e-4,si, w[id],&res[id],&err[id]);
				J(i,j) +=res[id];
				//cout << "  " << res << endl;
				//cout << GSL_FN_EVAL(&F,(grid(i)+grid(i+1))*0.55)<<endl;

			}
			

		}
	}
	//gsl_integration_workspace_free(w);
	for (int i=0; i < OMPTHREADS; i++){
		gsl_integration_workspace_free(w[i]);
		delete [] sing[i];
		delete [] param[i];
	}
	
}

