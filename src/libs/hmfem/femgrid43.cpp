#include "femgrid43.hpp"

using namespace HMFem::Impl;

Grid43::Grid43(GridGeom* parent){
	//Temp. implementation: only 4/3 grids on input
	assert([](GridGeom* g){
		for (int i=0; i<g->n_cells(); ++i){
			if (g->get_cell(i)->dim()>4) return false;
		}
		return true;
	}(parent));
	//shallow copy of points, deep copy of cells
	shallow_copy(parent, this);

	cells.clear();
	for (int i=0; i<parent->n_cells(); ++i){
		auto cell_old = parent->get_cell(i);
		auto c = aa::add_shared(cells, Cell());
		for (int j=0; j<cell_old->dim(); ++j){
			add_point_to_cell(c, const_cast<GridPoint*>(cell_old->get_point(i)));
		}
	}
}
