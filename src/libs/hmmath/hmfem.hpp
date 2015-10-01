#ifndef HYBMESH_HMFEM_HPP
#define HYBMESH_HMFEM_HPP

#include "hmproject.h"
#include "femgrid43.hpp"
#include "femmat.hpp"

namespace HMFem{

//Solution of (nabla^2)f = 0 equation.
class LaplasProblem{
	typedef std::function<double(const GridPoint*)> TDirFunc;
	typedef std::function<double(GridPoint*, GridPoint*)> TNeuFunc;
	struct TNeuData{
		GridPoint *point1, *point2;
		TNeuFunc *fun;
		double dist;
	};
	struct TDirData{
		const GridPoint* point;
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
	shared_ptr<HMFem::Grid43> grid;
	//Matricies
	shared_ptr<HMFem::Mat> laplas_mat;

	HMFem::Mat solution_mat;
	HMFem::MatSolve solver;
	vector<double> rhs;
	//boundary condition data
	std::list<TNeuFunc> _neufunc;
	std::list<TDirFunc> _dirfunc;
	std::set<TNeuData, TNeuCmp> neumann_data;
	std::set<TDirData, TDirCmp> dirichlet_data;

	//
	void RebuildSolutionMatrix();
public:
	//constructors:
	//creating grid and build laplas operator
	LaplasProblem(GridGeom* g);
	//using prebuild grid
	LaplasProblem(shared_ptr<Grid43> g);
	//using prebuilt grid and laplas matrix
	LaplasProblem(shared_ptr<Grid43> g, shared_ptr<Mat> lap);

	//boundary conditions: using vector of grid points
	void SetDirichlet(const vector<const GridPoint*>& pts, TDirFunc f);
	void SetNeumann(const vector<const GridPoint*>& pts, TNeuFunc f);

	//using contours which shear points with grid.
	void SetDirichlet(const HMCont2D::Contour& pts, TDirFunc f);

	//solve Ax=0 with bc rebuilding
	void Solve(vector<double>& ans);

	//solve Ax=0 without bc rebuilding. Could be called
	//only after Solve.
	void QuickSolve(vector<double>& ans);

	//f - is the solution of Laplas problem.
	//pnt should be sorted as an open path ( pnt[0]->pnt[back] )
	//Neumann conditions should exist on edges adjacent to this path.
	//Otherwise dfdn=0 on that edges is assumed.
	double IntegralDfDn(const vector<const GridPoint*>& pnt,
			const vector<double>& f);
	
	//pnt Contour should share points with grid
	double IntegralDfDn(const HMCont2D::Contour& pnt,
			const vector<double>& f);
};


}


#endif
