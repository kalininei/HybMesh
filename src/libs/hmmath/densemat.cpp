#include "densemat.hpp"
#include <numeric>

using namespace HMMath;

//lapack procedures: lu factorization and solution
extern "C" void dgetrf_(void*, void*, void*, void*, void*, void*);
extern "C" void dgetrs_(void*, void*, void*, void*, void*, void*, void*, void*, void*);

double DenseMat::RowMultVec(const vector<double>& u, int irow) const{
	return std::inner_product(dt.begin()+irow*N, dt.begin()+(irow+1)*N, u.begin(), 0.0);
}
void DenseMat::MultVec(const vector<double>& u, vector<double>& res) const{
	for (int i=0; i<N; ++i) res[i] += RowMultVec(u, i);
}

DenseSolver::DenseSolver(const DenseMat& mat){
	N = mat.N;
	lumat.resize(N*N);
	ipiv.resize(N);
	int INFO = 0;
	std::copy(mat.dt.begin(), mat.dt.end(), lumat.begin());
	dgetrf_(&N, &N, &lumat[0], &N, &ipiv[0], &INFO);
	if (INFO!=0) throw std::runtime_error("dense matrix solution initialization failed");
}

void DenseSolver::Solve(const vector<double>& rhs, vector<double>& x) const{
	int INFO = 0;
	char T = 'T'; //transpose before solution since lapack has fortran matrix ordering
	int Nrhs = 1;

	std::copy(rhs.begin(), rhs.begin()+N, x.begin());
	dgetrs_(&T, &N, &Nrhs, &lumat[0], &N, &ipiv[0], &x[0], &N, &INFO);
	if (INFO!=0) throw std::runtime_error("dense matrix solution failed");
}
