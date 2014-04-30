#ifndef DERIVATIVE_H
#define DERIVATIVE_H
#include "vecmat.h"
#include <vector>
#include <string>

class Derivative{

public:
	DEVector grid, tgrid;
	vector<DEVector> u;

	Derivative(DEVector grid_, DEVector tgrid_, vector<DEVector> u_)
	: grid(grid_), tgrid(tgrid_), u(u_){}
	
	virtual void initialize()=0;
	
	virtual void solve()=0;
	
	virtual void write_to_file(std::string filename)=0;
};

#endif
