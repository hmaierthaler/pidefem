#ifndef LEVY_FEM_H
#define LEVY_FEM_H

typedef double (*Func1D)(double);

void assemble_diff_op(T sigma, T gamma, T lambda, DEVector grid, GEMatrix &A);
void assemble_mass_mtrx(DEVector grid, GEMatrix &M);

void assemble_rhs_eur_put(	DEVector &rhs, DEVector grid, double offset,
		 	double tau, double sigma, double gamma, double lambda, double r, double K,
		 	double C, double G, double M, double Y);
		 	
void assemble_rhs_am_put(	DEVector &rhs, DEVector grid, double offset,
			double sigma, double gamma, double lambda, double r, double K,
			double C, double G, double M, double Y);

void assemble_jump_mat(DEVector grid, GEMatrix &J, double C,double G, double M, double Y);

void assemble_rhs_swing_am_put(	DEVector &rhs, DEVector grid, double offset,
				double tau, double delta, int p, 
				double sigma, double gamma, double lambda, double r, double K,
				Func1D C, double G, double M, double Y);
				
void assemble_rhs_swing_eur_put(	DEVector &rhs, DEVector grid, double offset,
				double tau, double delta, int p, 
				double sigma, double gamma, double lambda, double r, double K,
				Func1D C, double G, double M, double Y);
				
				
				
void assemble_rhs_eur_call(	DEVector &rhs, DEVector grid, double offset,
		 	double tau, double sigma, double gamma, double lambda, double r, double K,
		 	double C, double G, double M, double Y);
		 	
void assemble_rhs_am_call(	DEVector &rhs, DEVector grid, double offset,
			double sigma, double gamma, double lambda, double r, double K,
			double C, double G, double M, double Y);
			
void assemble_rhs_swing_am_call(	DEVector &rhs, DEVector grid, double offset,
				double tau, double delta, int p, 
				double sigma, double gamma, double lambda, double r, double K,
				Func1D C, double G, double M, double Y);
				
void assemble_rhs_swing_eur_call(	DEVector &rhs, DEVector grid, double offset,
				double tau, double delta, int p, 
				double sigma, double gamma, double lambda, double r, double K,
				Func1D C, double G, double M, double Y);
				


#endif
