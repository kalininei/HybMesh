#include "hmfem.hpp"
#include "femassembly.hpp"

using namespace HMFem;

LaplasProblem::LaplasProblem(GridGeom* g):
		grid(g),
		laplas_mat(HMFem::Impl::Assemble::PureLaplas(grid)),
		solution_mat(),
		rhs(grid.n_points(), 0.0){}

void LaplasProblem::SetDirichlet(const vector<GridPoint*>& pts, TDirFunc f){
	_THROW_NOT_IMP_;
}

void LaplasProblem::SetNeumann(const vector<GridPoint*>& pts, TNeuFunc f){
	_THROW_NOT_IMP_;
}

void LaplasProblem::RebuildSolutionMatrix(){
	std::fill(rhs.begin(), rhs.end(), 0.0);
	solution_mat = laplas_mat;

	//Neumann
	for (auto& nc: neumann_data){
		//val = (df/dn)*L/2;
		double val = ((*nc.fun)(nc.point1, nc.point2)) * nc.dist / 2;
		rhs[nc.point1->get_ind()] += val;
		rhs[nc.point2->get_ind()] += val;
	}

	//Dirichlet. Strictly after Neumann
	for (auto& dc: dirichlet_data){
		//put 1 to diagonal and value to rhs
		int ind = dc.point->get_ind();
		double val = (*dc.fun)(dc.point);
		solution_mat.clear_row(ind);
		solution_mat.set(ind, ind, 1.0);
		rhs[ind] = val;
	}

	//solver initialization
	solver.Init(solution_mat);
}


//solve Ax=0
void LaplasProblem::Solve(vector<double>& ans){
	RebuildSolutionMatrix();
	QuickSolve(ans);
}


void LaplasProblem::QuickSolve(vector<double>& ans){
	solver.Solve(rhs, ans);
}

