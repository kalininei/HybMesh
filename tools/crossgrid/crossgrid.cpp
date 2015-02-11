#include "crossgrid.h"
#include "grid.h"
#include "fileproc.h"

Grid* grid_construct(int Npts, int Ncells, double* pts, int* cells){
	return new GridGeom(Npts, Ncells, pts, cells);
}

int grid_npoints(Grid* g){
	return static_cast<GridGeom*>(g)->n_points();
}

int grid_ncells(Grid* g){
	return static_cast<GridGeom*>(g)->n_cells();
}
int grid_cellsdim(Grid* g){
	return static_cast<GridGeom*>(g)->n_cellsdim();
}

void grid_get_points_cells(Grid* g, double* pts, int* cells){
	auto gg=static_cast<GridGeom*>(g);
	//points
	for (int i=0; i<gg->n_points(); ++i){
		auto p = gg->get_point(i);
		*pts++ = p->x;
		*pts++ = p->y;
	}
	//cells
	for (int i=0; i<gg->n_cells(); ++i){
		auto c = gg->get_cell(i);
		for (int j=0; j<c->dim(); ++j){
			*cells++ = c->get_point(j)->get_ind();
		}
	}
}

void grid_save_vtk(Grid* g, const char* fn){
	save_vtk(static_cast<GridGeom*>(g), fn);
}

void grid_free(Grid* g){
	delete g;
}

Grid* cross_grids(Grid* gbase, Grid* gsecondary, double buffer_size){
	return GridGeom::cross_grids2(
			static_cast<GridGeom*>(gbase),
			static_cast<GridGeom*>(gsecondary),
			buffer_size);
}

