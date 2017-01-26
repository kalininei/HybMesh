#ifndef HYBMESH_HMFEM_HPP
#define HYBMESH_HMFEM_HPP

#include "hmproject.h"
#include "femgrid43.hpp"
#include "spmat.hpp"

namespace HMFem{

//Solution of (nabla^2)f = 0 equation.
class LaplasProblem{
	typedef std::function<double(const HM2D::Vertex*)> TDirFunc;
	typedef std::function<double(const HM2D::Vertex*, const HM2D::Vertex*)> TNeuFunc;
	struct TNeuData{
		int index1, index2;
		TNeuFunc *fun;
		double dist;
	};
	struct TDirData{
		int index;
		TDirFunc* fun;
	};
	struct TNeuCmp{
		bool operator()(const TNeuData& x, const TNeuData& y) const {
			if (x.index1 == y.index1){
				return x.index2 < y.index2;
			} else { 
				return x.index1 < y.index1;
			}
		}
	};
	struct TDirCmp{
		bool operator()(const TDirData& x, const TDirData& y) const {
			return x.index < y.index;
		}
	};

	//Grids
	const HM2D::GridData* grid;
	//Matricies
	shared_ptr<HMMath::Mat> laplas_mat;

	HMMath::Mat solution_mat;
	shared_ptr<HMMath::MatSolve> solver;
	vector<double> rhs;
	//boundary condition data
	std::list<TNeuFunc> _neufunc;
	std::list<TDirFunc> _dirfunc;
	std::set<TNeuData, TNeuCmp> neumann_data;
	std::set<TDirData, TDirCmp> dirichlet_data;

	void RebuildSolutionMatrix();

	mutable std::set<HM2D::Vertex> _bp;
	const HM2D::Vertex* get_boundary_point(const Point& p) const;
	int get_boundary_point_index(const Point& p) const;
public:
	//constructors:
	//build laplas operator
	LaplasProblem(const HM2D::GridData& g);
	//using prebuilt laplas matrix
	LaplasProblem(const HM2D::GridData& g, shared_ptr<HMMath::Mat> lap);

	//boundary conditions: using vector of grid points
	void ClearBC();
	void SetDirichlet(const HM2D::VertexData& pts, TDirFunc f);
	void SetDirichlet(const vector<int>& pts, TDirFunc f);
	void SetNeumann(const HM2D::VertexData& pts, TNeuFunc f);

	//using contours which shear points with grid.
	void SetDirichlet(const HM2D::EdgeData& pts, TDirFunc f);

	//solve Ax=0 with bc rebuilding
	void Solve(vector<double>& ans);

	//solve Ax=0 without bc rebuilding. Could be called
	//only after Solve.
	void QuickSolve(vector<double>& ans);

	//solve Ax=0 with rebuilding of rhs vector with bc
	void QuickSolve_BC(vector<double>& ans);

	//f - is the solution of Laplas problem.
	//pnt should be sorted as an open path ( pnt[0]->pnt[back] )
	//Neumann conditions should exist on edges adjacent to this path.
	//Otherwise dfdn=0 on that edges is assumed.
	double IntegralDfDn(const vector<const HM2D::Vertex*>& pnt,
			const vector<double>& f);
	
	//pnt Contour should share points with grid
	double IntegralDfDn(const HM2D::EdgeData& pnt,
			const vector<double>& f);
};


}


#endif
