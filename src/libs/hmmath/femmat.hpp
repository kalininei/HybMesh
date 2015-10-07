#ifndef HYBMESH_HMFEM_FEMMAT_HPP
#define HYBMESH_HMFEM_FEMMAT_HPP
#include "hmproject.h"
#include "femgrid43.hpp"

namespace HMFem{

//Sparce Matrix
struct Mat{
	vector<std::map<int, double>> data;

	//get/set methods
	int rows() const {return data.size();}
	void clear_row(int i) { data[i].clear(); }
	void set(int i, int j, double val);
	double get(int i, int j) const;
	void add(int i, int j, double val);
	vector<double> diag() const;
	int row_size(int irow) const { return data[irow].size(); }
	int nnz() const {
		return std::accumulate(data.begin(), data.end(), 0,
			[](int a, const std::map<int, double>& m){
				return a+m.size();
			});
	}

	//Methods
	double RowMultVec(const vector<double>& u, int irow) const;
	void MultVec(const vector<double>& u, vector<double>& res) const;
};


//dense matrix of little dimension
struct LocMat{
	virtual void ToMat(const vector<int>& pind, Mat& target) const = 0;
};

struct LocMat3Sym: public LocMat, public std::array<double, 6>{
	void ToMat(const vector<int>& pind, Mat& target) const;
};

struct LocMat3: public LocMat, public std::array<double, 9>{
	void ToMat(const vector<int>& pind, Mat& target) const;
};

struct LocMat4Sym: public LocMat, public std::array<double, 10>{
	void ToMat(const vector<int>& pind, Mat& target) const;
};

// ==================== Solving procedures
class MatSolve{
protected:
	const Mat* m;
public:
	MatSolve(const Mat& mat): m(&mat){}
	virtual ~MatSolve(){}
	struct Options{
		Options():
			direct_solver_max_nnz(200000),
			iter_tol(geps),
			iter_maxit(10000)
		{}
		int direct_solver_max_nnz;
		double iter_tol;
		int iter_maxit;
	};
	virtual void Solve(const vector<double>& rhs, vector<double>& x) = 0;

	static shared_ptr<MatSolve>
	Factory(const Mat& m, Options=Options());
};

class SeidelSolver: public MatSolve{
	int MaxIt;
	double tol;
public:
	SeidelSolver(const Mat& m, Options opt);
	void Solve(const vector<double>& rhs, vector<double>& x) override;
};

class SuiteSparseQRSolver: public MatSolve{
	void* slv;
public:
	SuiteSparseQRSolver(const Mat& m, Options opt);
	~SuiteSparseQRSolver();
	void Solve(const vector<double>& rhs, vector<double>& x) override;
};



}
#endif
