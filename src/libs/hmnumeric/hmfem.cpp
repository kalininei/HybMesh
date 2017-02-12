#include "hmfem.hpp"
#include "femassembly.hpp"
#include "assemble2d.hpp"

using namespace HMFem;

const HM2D::Vertex* LaplasProblem::get_boundary_point(const Point& p) const{
	int ind = get_boundary_point_index(p);
	return grid->vvert[ind].get();
}
int LaplasProblem::get_boundary_point_index(const Point& p) const{
	if (_bp.size() == 0){
		auto be = HM2D::ECol::Assembler::GridBoundary(*grid);
		auto av = HM2D::AllVertices(be);
		aa::enumerate_ids_pvec(grid->vvert);
		for (auto pbp: av){
			_bp.insert(*pbp);
		}
	}
	auto fnd = _bp.find(p);
	return fnd->id;
}

LaplasProblem::LaplasProblem(const HM2D::GridData& g):
		grid(&g),
		laplas_mat(HMFem::Assemble::PureLaplas(*grid)),
		solution_mat(),
		rhs(grid->vvert.size(), 0.0){}

LaplasProblem::LaplasProblem(const HM2D::GridData& g, shared_ptr<HMMath::Mat> lap):
		grid(&g), laplas_mat(lap),
		solution_mat(),
		rhs(grid->vvert.size(), 0.0){}

void LaplasProblem::ClearBC(){
	_neufunc.clear();
	_dirfunc.clear();
	neumann_data.clear();
	dirichlet_data.clear();
}
void LaplasProblem::SetDirichlet(const HM2D::VertexData& pts, TDirFunc f){
	_dirfunc.push_back(f);
	for (auto p: pts){
		int ind = get_boundary_point_index(*p);
		TDirData dt {ind, &_dirfunc.back()};
		dirichlet_data.insert(dt);
	}
}
void LaplasProblem::SetDirichlet(const vector<int>& pts, TDirFunc f){
	_dirfunc.push_back(f);
	for (auto ind: pts){
		TDirData dt {ind, &_dirfunc.back()};
		dirichlet_data.insert(dt);
	}
}
void LaplasProblem::SetDirichlet(const HM2D::EdgeData& pts, TDirFunc f){
	SetDirichlet(HM2D::AllVertices(pts), f);
}

void LaplasProblem::SetNeumann(const HM2D::VertexData& pts, TNeuFunc f){
	_THROW_NOT_IMP_;
}

void LaplasProblem::RebuildSolutionMatrix(){
	std::fill(rhs.begin(), rhs.end(), 0.0);
	solution_mat = *laplas_mat;

	//Neumann
	for (auto& nc: neumann_data){
		//val = (df/dn)*L/2;
		double val = ((*nc.fun)(grid->vvert[nc.index1].get(), grid->vvert[nc.index2].get())) * nc.dist / 2;
		rhs[nc.index1] += val;
		rhs[nc.index2] += val;
	}

	//Dirichlet. Strictly after Neumann
	for (auto& dc: dirichlet_data){
		//put 1 to diagonal and value to rhs
		double val = (*dc.fun)(grid->vvert[dc.index].get());
		solution_mat.clear_row(dc.index);
		solution_mat.set(dc.index, dc.index, 1.0);
		rhs[dc.index] = val;
	}

	//solver initialization
	solver = HMMath::MatSolve::Factory(solution_mat);
}

//solve Ax=0
void LaplasProblem::Solve(vector<double>& ans){
	RebuildSolutionMatrix();
	QuickSolve(ans);
}


void LaplasProblem::QuickSolve(vector<double>& ans){
	solver->Solve(rhs, ans);
}

void LaplasProblem::QuickSolve_BC(vector<double>& ans){
	std::fill(rhs.begin(), rhs.end(), 0.0);

	//Neumann
	for (auto& nc: neumann_data){
		//val = (df/dn)*L/2;
		double val = ((*nc.fun)(grid->vvert[nc.index1].get(), grid->vvert[nc.index2].get())) * nc.dist / 2;
		rhs[nc.index1] += val;
		rhs[nc.index2] += val;
	}

	//Dirichlet. Strictly after Neumann
	for (auto& dc: dirichlet_data){
		//put 1 to diagonal and value to rhs
		double val = (*dc.fun)(grid->vvert[dc.index].get());
		solution_mat.clear_row(dc.index);
		solution_mat.set(dc.index, dc.index, 1.0);
		rhs[dc.index] = val;
	}

	solver->Solve(rhs, ans);
}

double LaplasProblem::IntegralDfDn(const vector<const HM2D::Vertex*>& pnt,
		const vector<double>& f){
	//we assume that pnt are ordered in such a way that
	//each pair <pnt[i], pnt[i+1]> forms an edge
	assert(pnt.size()>1 && pnt[0] != pnt.back());
	assert( [&](){
		auto tree = HM2D::Contour::Tree::GridBoundary(*grid);
		shared_ptr<HM2D::Contour::Tree::TNode> c;
		for (auto& t: tree.nodes){
			for (auto v: HM2D::AllVertices(t->contour)){
				if (*pnt[0] == *v){
					c = t;
					break;
				}
			}
			if (c) break;
		}
		if (!c) return false;
		auto pp = HM2D::Contour::OrderedPoints(c->contour); pp.pop_back();
		auto fnd = pp.begin();
		while (**fnd != *pnt[0]) ++fnd;
		//choose forward or backward
		std::rotate(pp.begin(), fnd, pp.end());
		if (*pnt[1] != *pp[1]){
		//try backward
			std::reverse(pp.begin(), pp.end());
			std::rotate(pp.begin(), pp.end()-1, pp.end());
		}
		for (int i=0; i<pnt.size(); ++i){
			if (*pnt[i] != *pp[i]) return false;
		}
		return true;
	}());
	aa::enumerate_ids_pvec(grid->vvert);
	
	//assemble pure rhs
	vector<double> r(pnt.size(), 0.0);
	for (int i=0; i<pnt.size(); ++i){
		const HM2D::Vertex* p1 = pnt[i];
		r[i] = laplas_mat->RowMultVec(f, p1->id);
	}
	//add neumann condition to first and last rhs entries
	//if it exists
	if (neumann_data.size()>0){
		_THROW_NOT_IMP_;
	}
	
	//find |L|*k*dfdn using simple elimination
	vector<double> lk_dfdn(pnt.size()-1, 0.0);
	lk_dfdn[0] = 2*r[0];
	for (int i=1; i<pnt.size()-1; ++i){
		lk_dfdn[i] = 2*r[i] - lk_dfdn[i-1];
	}

	//return sum
	return std::accumulate(lk_dfdn.begin(), lk_dfdn.end(), 0.0);
}

double LaplasProblem::IntegralDfDn(const HM2D::EdgeData& pnt, const vector<double>& f){
	vector<const HM2D::Vertex*> pts2; pts2.reserve(pnt.size());
	for(auto p: HM2D::Contour::OrderedPoints(pnt))
		pts2.push_back(get_boundary_point(*p));
	return IntegralDfDn(pts2, f);
}
