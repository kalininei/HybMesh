#ifndef HYBMESH_HMFEM_HPP
#define HYBMESH_HMFEM_HPP

#include "hmproject.h"
#include "femgrid43.hpp"
#include "femmat.hpp"

namespace HMFem{

class LaplasProblem{
	typedef std::function<double(GridPoint*)> TDirFunc;
	typedef std::function<double(GridPoint*, GridPoint*)> TNeuFunc;
	struct TNeuData{
		GridPoint *point1, *point2;
		TNeuFunc *fun;
		double dist;
	};
	struct TDirData{
		GridPoint* point;
		TDirFunc* fun;
	};
	struct TNeuCmp{
		bool operator()(const TNeuData& x, const TNeuData& y) const {
			return (x.point1->get_ind() == y.point1->get_ind()) ?
				x.point2->get_ind() < y.point2->get_ind() :
				x.point1->get_ind() < y.point1->get_ind();
		}
	};
	struct TDirCmp{
		bool operator()(const TDirData& x, const TDirData& y) const {
			return x.point->get_ind() < y.point->get_ind();
		}
	};

	//Grids
	HMFem::Impl::Grid43 grid;
	//Matricies
	HMFem::Impl::Mat laplas_mat;
	HMFem::Impl::Mat solution_mat;
	HMFem::Impl::MatSolve solver;
	vector<double> rhs;
	//boundary condition data
	std::set<TNeuFunc> _neufunc;
	std::set<TDirFunc> _dirfunc;
	std::set<TNeuData, TNeuCmp> neumann_data;
	std::set<TDirData, TDirCmp> dirichlet_data;

	//
	void RebuildSolutionMatrix();
public:
	//constructor
	LaplasProblem(GridGeom* g);

	//boundary conditions
	void SetDirichlet(const vector<GridPoint*>& pts, TDirFunc f);
	void SetNeumann(const vector<GridPoint*>& pts, TNeuFunc f);

	//solve Ax=0 with bc rebuilding
	void Solve(vector<double>& ans);

	//solve Ax=0 without bc rebuilding
	void QuickSolve(vector<double>& ans);
	
};


}


#endif
