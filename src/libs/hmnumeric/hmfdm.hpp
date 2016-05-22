#ifndef HMMATH_HMFDM_HPP
#define HMMATH_HMFDM_HPP
#include "hmproject.h"
#include "spmat.hpp"

namespace HMFdm{

//laplas solver is square area with regular mesh
class LaplasSolver{
	vector<double> x, y;
	vector<double> rhs;
	std::map<int, double> predefined_values;
	void set_predef_value(int i, int j, double val);
	shared_ptr<HMMath::MatSolve> solver;
	bool was_init(){ return solver != nullptr; }
	void initialize(); 
	void assemble_rhs();
public:
	// ==== properties
	int N() const { return x.size()*y.size(); }
	int Nx() const { return x.size(); }
	int Ny() const { return y.size(); }
	int glob_index(int i, int j) const { return j*x.size() + i; }
	std::pair<int, int> sub_index(int gi) const { return std::make_pair<int, int>(gi%Nx(), gi/Nx()); }

	// ==== conststructor
	LaplasSolver(const vector<double>& _x, const vector<double>& _y): x(_x), y(_y){}

	// ==== set boundary conditions
	enum class Bnd {All, Top, Bottom, Left, Right};
	void SetBndValues(Bnd b, const std::function<double(int, int)>& f);

	// ==== Solver
	//solves problem: Laplas(ans) = 0.
	//could be called multipole times with changed boundary values,
	//but boundary stencil could not be changed after first Solve call.
	void Solve(vector<double>& ans);
};

};

#endif
