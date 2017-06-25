#ifndef HYBMESH_DENSEMAT_HPP
#define HYBMESH_DENSEMAT_HPP

#include "hmproject.h"

namespace HMMath{

//NxN matrix
struct DenseMat{
	const int N;
	std::vector<double> dt;

	DenseMat(int N): N(N), dt(N*N){};

	double val(int i, int j) const { return dt[i*N+j]; }
	double& val(int i, int j) { return dt[i*N+j]; }

	//Methods
	double RowMultVec(const vector<double>& u, int irow) const;
	void MultVec(const vector<double>& u, vector<double>& res) const;
};


class DenseSolver{
	mutable vector<int> ipiv;
	mutable vector<double> lumat;
	mutable int N;
public:
	DenseSolver(const DenseMat& mat);
	void Solve(const vector<double>& rhs, vector<double>& x) const;
};

}
#endif
