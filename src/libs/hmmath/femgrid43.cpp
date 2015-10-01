#include "femgrid43.hpp"
#include "trigrid.h"

using namespace HMFem;

shared_ptr<Grid43>
Grid43::Build(GridGeom* parent){
	shared_ptr<Grid43> ret(new Grid43());
	//Temp. implementation: only 4/3 grids on input
	assert([](GridGeom* g){
		for (int i=0; i<g->n_cells(); ++i){
			if (g->get_cell(i)->dim()>4) return false;
		}
		return true;
	}(parent));
	//shallow copy of points, deep copy of cells
	shallow_copy(parent, ret.get());

	ret->cells.clear();
	for (int i=0; i<parent->n_cells(); ++i){
		auto cell_old = parent->get_cell(i);
		auto c = aa::add_shared(ret->cells, Cell());
		for (int j=0; j<cell_old->dim(); ++j){
			add_point_to_cell(c, const_cast<GridPoint*>(cell_old->get_point(i)));
		}
	}

	return ret;
}

shared_ptr<Grid43>
Grid43::Build3(const vector<Point>& p, double h){
	shared_ptr<Grid43> ret(new Grid43());
	shared_ptr<TriGrid> g1 = TriGrid::TriangulateArea(p, h);
	GGeom::Modify::ShallowAdd(g1.get(), ret.get());
	return ret;
}

shared_ptr<Grid43>
Grid43::Build3(GridGeom* orig){
	shared_ptr<Grid43> ret(new Grid43());
	//copy points and cells
	GGeom::Modify::ShallowAdd(orig, ret.get());
	//remove cells
	ret->cells.clear();
	//build triangle cells
	auto addtri = [&](Cell* oldcell, int i1, int i2, int i3){
		Cell* newcell = aa::add_shared(ret->cells, Cell());
		add_point_to_cell(newcell, const_cast<GridPoint*>(oldcell->get_point(i1)));
		add_point_to_cell(newcell, const_cast<GridPoint*>(oldcell->get_point(i2)));
		add_point_to_cell(newcell, const_cast<GridPoint*>(oldcell->get_point(i3)));
	};
	for (int i=0; i<orig->n_cells(); ++i){
		auto oc = orig->get_cell(i);
		if (oc->dim() == 3) addtri(oc, 0, 1, 2);
		else if (oc->dim() == 4){ addtri(oc, 0, 1, 2); addtri(oc, 0, 2, 3); }
		else _THROW_NOT_IMP_;
	}
	//set indicies
	ret->set_indicies();
	return ret;
}

shared_ptr<Grid43::Approximator> Grid43::GetApprox() const{
	return shared_ptr<Grid43::Approximator>(new Grid43::Approximator(this));
}


double Grid43::Approximator::Val(Point p, const vector<double>& fun) const{
	_THROW_NOT_IMP_;
}


vector<double> Grid43::Approximator::Vals(Point, const vector<const vector<double>*>& funs) const{
	_THROW_NOT_IMP_;
}
