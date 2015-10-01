#include "femmat.hpp"

using namespace HMFem;

namespace{

void SeidelSolve(const Mat& M, const vector<double>& rhs, vector<double>& u){
	const double tol = 1e-9;
	const int maxit = 10000;

	vector<double> diag = M.diag();
	assert(std::all_of(diag.begin(), diag.end(), [](double x){ return x>0; }));
	for (int it=0; it<maxit; ++it){
		double rNorm = 0, r=0;
		for (int i=0; i<M.rows(); ++i){
			r = rhs[i] - M.RowMultVec(u, i);
			u[i] += r/diag[i];
			if (fabs(r)>rNorm) rNorm = fabs(r);
		}
		if (rNorm<=tol) return;
	}
	throw std::runtime_error("Seidel solver failed to converge");
}

}

void Mat::set(int i, int j, double val) {
	if (i >= data.size()) data.resize(i+1);
	data[i][j] = val;
}

void Mat::add(int i, int j, double val) {
	if (i >= data.size()) data.resize(i+1);
	auto fnd = data[i].emplace(j, 0.0);
	fnd.first->second += val;
}

double Mat::get(int i, int j) const{
	if (i >= data.size()) return 0;
	auto fnd = data[i].find(j);
	if (fnd == data[i].end()) return 0;
	else return fnd->second;
};

vector<double> Mat::diag() const {
	vector<double> ret(rows());
	int i=0;
	std::transform(data.begin(), data.end(), ret.begin(), [&i](const std::map<int, double>& dt){
		auto fnd = dt.find(i++);
		if (fnd == dt.end()) return 0.0;
		else return fnd->second;
	});
	return ret;
}

double Mat::RowMultVec(const vector<double>& u, int irow) const {
	assert(irow<data.size());

	return std::accumulate(data[irow].begin(), data[irow].end(), 0.0,
		[&u](double s, const std::pair<int, double>& it){
			return s + it.second*u[it.first];
		}
	);
}


void MatSolve::Init(Mat& mat){
	m = &mat;
}

void MatSolve::Solve(const vector<double>& rhs, vector<double>& x){
	//Temporary using Seidel solver
	SeidelSolve(*m, rhs, x);
}


void LocMat3Sym::ToMat(const vector<int>& pind, Mat& target) const {
	target.add(pind[0], pind[0], (*this)[0]);
	target.add(pind[0], pind[1], (*this)[1]);
	target.add(pind[0], pind[2], (*this)[2]);
	target.add(pind[1], pind[0], (*this)[1]);
	target.add(pind[2], pind[0], (*this)[2]);

	target.add(pind[1], pind[1], (*this)[3]);
	target.add(pind[1], pind[2], (*this)[4]);
	target.add(pind[2], pind[1], (*this)[4]);

	target.add(pind[2], pind[2], (*this)[5]);
}

void LocMat4Sym::ToMat(const vector<int>& pind, Mat& target) const {
	target.add(pind[0], pind[0], (*this)[0]);
	target.add(pind[0], pind[1], (*this)[1]);
	target.add(pind[0], pind[2], (*this)[2]);
	target.add(pind[0], pind[3], (*this)[3]);
	target.add(pind[1], pind[0], (*this)[1]);
	target.add(pind[2], pind[0], (*this)[2]);
	target.add(pind[3], pind[0], (*this)[3]);

	target.add(pind[1], pind[1], (*this)[4]);
	target.add(pind[1], pind[2], (*this)[5]);
	target.add(pind[1], pind[3], (*this)[6]);
	target.add(pind[2], pind[1], (*this)[5]);
	target.add(pind[3], pind[1], (*this)[6]);

	target.add(pind[2], pind[2], (*this)[7]);
	target.add(pind[2], pind[3], (*this)[8]);
	target.add(pind[3], pind[2], (*this)[8]);

	target.add(pind[3], pind[3], (*this)[9]);
}



