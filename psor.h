#ifndef PSOR_H
#define PSOR_H 1

#include <cmath>


template <typename MATTYPE, typename VECTYPE>
    int
    psor(MATTYPE &A, VECTYPE &x,VECTYPE x0, const VECTYPE &b, const VECTYPE &c, double tol = 0.0000000000001, int maxIterations = 20000);


    //sample call: psor(A,x,b,c);

    template <typename MATTYPE, typename VECTYPE>
    int
    psor(MATTYPE &A, VECTYPE &x, VECTYPE x0, const VECTYPE &b, const VECTYPE &c, double tol, int maxIterations){
    	double omega=0.95,
    		v,w;
    	int n = A.numCols();
        VECTYPE xp=x0,dx;
	

        if ((int)x.length()!=A.numCols()) {
            x = VECTYPE(A.numCols());
        }

	
        for (int k=1; k<=maxIterations; k++) {
        
        	for (int i=1; i<=n; i++){
        		v=0;
        		w=0;
        	
        		for (int j=1; j <= i-1; j++){
        			v += A(i,j)*x(j);
        		}
        		for (int j=i+1; j<= n; j++){
        			w += A(i,j)*xp(j);
        		}
        		x(i) = (b(i) - v - w)/A(i,i);
        		x(i) = std::max(c(i), xp(i) + omega*(x(i)-xp(i))) ;
        		//x(i)=xp(i) + omega*(x(i)-xp(i));


        		//xp(i) =x(i);
        	}
        	dx=x-xp;
        	if ( dx*dx < tol*tol){
        		return k;
        	}
        	xp=x;
        }

        return maxIterations;
    }


#endif // PSOR_H

