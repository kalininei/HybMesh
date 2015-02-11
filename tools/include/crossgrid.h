#ifndef CROSSGRID_CROSSGRID_H
#define CROSSGRID_CROSSGRID_H

//base grid class
class Grid{
public:
	virtual ~Grid(){}
};

extern "C"{

	Grid* grid_construct(int Npts, int Ncells, double* pts, int* cells);

	int grid_npoints(Grid* g);
	int grid_ncells(Grid* g);
	int grid_cellsdim(Grid* g);
	void grid_get_points_cells(Grid* g, double* pts, int* cells);
	void grid_save_vtk(Grid* g, const char* fn);

	void grid_free(Grid* g);

	Grid* cross_grids(Grid* gbase, Grid* gsecondary, double buffer_size);

};


#endif
