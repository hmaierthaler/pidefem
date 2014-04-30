#include "vecmat.h"

// assemble diffusion part of the stiffness matrix
void assemble_diff_op(T sigma, T gamma, T lambda, DEVector grid, GEMatrix &A){

	A=0;
	for (int i=1; i <= grid.length(); i++){
	
		
		if ( i > 1){


			A(i,i) = -0.5*sigma*sigma /(grid(i)-grid(i-1));
			A(i,i-1) = -0.5*sigma*sigma /(grid(i-1)-grid(i)) 
				- 0.5*(gamma) ;
			//OU Part
			A(i,i) += -lambda*(grid(i-1)/6.0+grid(i)/3.0);
			A(i,i-1) += lambda*(grid(i)/3.0 +grid(i-1)/6.0);
		}
		
		if ( i < grid.length() ){

			
			A(i,i) += -0.5*sigma*sigma /(grid(i+1)-grid(i));
			A(i,i+1) += -0.5*sigma*sigma /(grid(i)-grid(i+1))
				 +0.5* (gamma) ;
			//OU Part
			A(i,i) += lambda*(grid(i)/3.0+grid(i+1)/6.0);
			A(i,i+1) += -lambda*(grid(i)/3.0 + grid(i+1)/6.0);
		}

	}
	
	/*
	// infinite element correction
	double a=2; //a=1/L
	A(1,1) += -0.5*sigma*sigma*0.5*a + 0.5*gamma;
	A(grid.length(),grid.length()) += -0.5*sigma*sigma*0.5*a -0.5*gamma;
	// OU
	A(1,1) +=-lambda*(-0.25/a +0.5*grid(1));
	A(grid.length(),grid.length()) += -lambda*(0.25/a -0.5*grid(grid.length()));
	*/
}

// assemble mass matrix
void assemble_mass_mtrx(DEVector grid, GEMatrix &M){

	for (int i=1; i <= grid.length(); i++){
		if ( i > 1){
			M(i,i)=(grid(i)-grid(i-1))/3.0;
			M(i,i-1) = (grid(i) - grid(i-1))/6.0;
		}
		
		if ( i < grid.length() ){
			M(i,i) += (grid(i+1) - grid(i))/3.0;
			M(i,i+1) += (grid(i+1) - grid(i))/6.0;

		}
	}
	
	/*
	//infinite element correction
	double a=2;
	M(1,1) += 0.5/a;
	M(grid.length(),grid.length()) += 0.5/a;
	*/
	

}
