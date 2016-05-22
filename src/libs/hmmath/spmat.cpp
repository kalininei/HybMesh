#include "spmat.hpp"
#include "cholmod.h"
#include "umfpack.h"
#include "SuiteSparseQR_C.h"

using namespace HMMath;

std::ostream& HMMath::operator<<(std::ostream& os, const Mat& mat){
	os<<"N = "<<mat.data.size()<<"; Nnz = "<<mat.nnz()<<std::endl;
	int i = 0;
	for (auto& d: mat.data){
		os<<"row = "<<i<<std::endl;
		for (auto kv: d){
			os<<"\t"<<kv.first<<" -> "<<kv.second<<std::endl;
		}
		++i;
	}
	return os;
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

void Mat::MultVec(const vector<double>& u, vector<double>& res) const{
	for (int i=0; i<data.size(); ++i) res[i] = RowMultVec(u, i);
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

void LocMat3::ToMat(const vector<int>& pind, Mat& target) const{
	target.add(pind[0], pind[0], (*this)[0]);
	target.add(pind[0], pind[1], (*this)[1]);
	target.add(pind[0], pind[2], (*this)[2]);
	target.add(pind[1], pind[0], (*this)[3]);
	target.add(pind[1], pind[1], (*this)[4]);
	target.add(pind[1], pind[2], (*this)[5]);
	target.add(pind[2], pind[0], (*this)[6]);
	target.add(pind[2], pind[1], (*this)[7]);
	target.add(pind[2], pind[2], (*this)[8]);
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


// ====================== Matrix Solvers
shared_ptr<MatSolve>
MatSolve::Factory(const Mat& m, Options opt){
	shared_ptr<MatSolve> ret;
	if (m.nnz() < opt.direct_solver_max_nnz){
		ret.reset(new SuiteSparseQRSolver(m, opt));
	} else {
		//ret.reset(new SeidelSolver(m, opt));
		_THROW_NOT_IMP_;
	}
	return ret;
}

// ======================= Seidel solver
SeidelSolver::SeidelSolver(const Mat& m, Options opt): MatSolve(m){
	MaxIt = opt.iter_maxit;
	tol = opt.iter_tol;
}

void SeidelSolver::Solve(const vector<double>& rhs, vector<double>& u){
	vector<double> diag = m->diag();
	assert(std::all_of(diag.begin(), diag.end(), [](double x){ return x>0; }));
	u.resize(diag.size(), 0.0);
	for (int it=0; it<MaxIt; ++it){
		double rNorm = 0, r=0;
		for (int i=0; i<m->rows(); ++i){
			r = rhs[i] - m->RowMultVec(u, i);
			u[i] += r/diag[i];
			if (fabs(r)>rNorm) rNorm = fabs(r);
		}
		if (rNorm<=tol) return;
	}
	throw std::runtime_error("Seidel solver failed to converge");
}

// ======================= SuiteSparse QR Solver
namespace QRImpl{

struct SPQR{
	SPQR(): A(0), x(0), b(0), QR(0){ cholmod_l_start(&common); }
	~SPQR(){
		if (A!=0) {cholmod_l_free_sparse(&A, &common); A=0;}
		if (x!=0) {cholmod_l_free_dense(&x, &common); x=0;}
		if (b!=0) {cholmod_l_free_dense(&b, &common); b=0;}
		if (y!=0) {cholmod_l_free_dense(&y, &common); y=0;}
		if (QR!=0) {SuiteSparseQR_C_free(&QR, &common); QR=0;}
		cholmod_l_finish(&common);
	}
	void InitMat(const Mat& m){
		cholmod_triplet *T = cholmod_l_allocate_triplet(
			m.rows(), m.rows(), m.nnz(), 0, CHOLMOD_REAL, &common);
		int irow = 0, cur = 0;
		for (auto& row: m.data){
			int nrow = row.size();
			//rows
			std::fill_n((SuiteSparse_long*)T->i + cur,
					nrow, irow);
			//columns and values
			for(auto& cv: row){
				((SuiteSparse_long*)T->j)[cur] = cv.first;
				((double*)T->x)[cur] = cv.second;
				++cur;
			}
			++irow;
		}
		T->nnz = cur;
		A = cholmod_l_triplet_to_sparse(T, 0, &common);
		cholmod_l_free_triplet(&T, &common);
	}
	void InitSlv(){
		b = cholmod_l_zeros(A->nrow, 1, A->xtype, &common);
		QR = SuiteSparseQR_C_factorize(SPQR_ORDERING_DEFAULT, SPQR_DEFAULT_TOL, A, &common);
		x = y = 0;  //will be allocated within solve
	}
	int Solve(const double* rhs, double* v){
		std::copy(rhs,rhs+A->nrow,(double*)b->x);
	
		// Y = Q'*B
		y = SuiteSparseQR_C_qmult(SPQR_QTX, QR, b, &common) ;
		if (y == 0) return 0;

		// X = R\(E*Y)
		x = SuiteSparseQR_C_solve(SPQR_RETX_EQUALS_B, QR, y, &common) ;
		if (x == 0) return 0;

		//copy solution
		std::copy((double*)x->x, (double*)x->x+A->nrow, v);
		cholmod_l_free_dense(&y, &common); y = 0;
		cholmod_l_free_dense(&x, &common); x = 0;
		
		return 1;
	}

	cholmod_common common;
	cholmod_sparse *A;
	cholmod_dense *x, *b, *y;
	SuiteSparseQR_C_factorization *QR;
};

}

SuiteSparseQRSolver::SuiteSparseQRSolver(const Mat& m, Options opt): MatSolve(m){
	auto slv1 = new QRImpl::SPQR();
	slv1->InitMat(m);
	slv1->InitSlv();
	slv = slv1;
}
SuiteSparseQRSolver::~SuiteSparseQRSolver(){
	delete static_cast<QRImpl::SPQR*>(slv);
}

void SuiteSparseQRSolver::Solve(const vector<double>& rhs, vector<double>& x){
	int ans = static_cast<QRImpl::SPQR*>(slv)->Solve(&rhs[0], &x[0]);
	if (ans != 1) throw std::runtime_error("QRSolver matrix solution failed");
}


