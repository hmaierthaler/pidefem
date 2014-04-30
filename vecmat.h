#ifndef VECMAT_H
#define VECMAT_H

#include <flens/flens.cxx>

using namespace std;
using namespace flens;


// definitions of matrix and vector types
typedef double					T;
typedef DenseVector<Array<T> >			DEVector;
typedef GeMatrix<FullStorage<T, ColMajor > >	GEMatrix;
typedef DenseVector<Array<int> >		IndexVector;

#endif
