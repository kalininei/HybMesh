#ifndef CROSSGRID_CROSSGRID_H
#define CROSSGRID_CROSSGRID_H

//base grid class
class Grid{
public:
	virtual ~Grid(){}
};

extern "C"{

	// === grid management
	Grid* grid_construct(int Npts, int Ncells, double* pts, int* cells);
	//grid structures count
	int grid_npoints(Grid* g);
	int grid_ncells(Grid* g);
	int grid_cellsdim(Grid* g);
	//get grid in points, edges->points, cells->edges format
	void grid_get_edges_info(Grid* g, int* Npnt, int* Neds, int* Ncls,
			double** pts,
			int** ed_pt,
			int** cls_dims,
			int** cls_eds);
	//free edges data
	void grid_free_edges_info(double** pts, int** ed_pt, int** cls_dims, int** cls_eds);
	//get grid in points, cells->points format. arrays should be allocated already
	void grid_get_points_cells(Grid* g, double* pts, int* cells);
	//save grid to vtk
	void grid_save_vtk(Grid* g, const char* fn);
	//free grid
	void grid_free(Grid* g);

	// === set default callbacks
	//crossgrid callback pointer is the function of the following arguments:
	//	(global procedure name, local procedure name, global percentage, local percentage)
	//percentages are doubles in [0, 1] range. If it is less then, then percentages of procedures is not tracking.
	//normally returns CALLBACK_OK. Should return CALLBACK_CANCEL for cancellation require
	const int CALLBACK_OK = 0;
	const int CALLBACK_CANCEL = 1;
	typedef int (*crossgrid_callback)(const char*, const char*, double, double);
	void crossgrid_cout_callback();  //callback to std::cout
	void crossgrid_silent_callback(); //no callback at all
	void crossgrid_set_callback(crossgrid_callback fun); //user defined callback

	// === crossgrid procedure
	//with callback defined globally
	Grid* cross_grids(Grid* gbase, Grid* gsecondary, double buffer_size, 
			double density, int preserve_bp);
	//with specified callback function
	Grid* cross_grids_wcb(Grid* gbase, Grid* gsecondary, double buffer_size, 
			double density, int preserve_bp, crossgrid_callback cb_fun);

	// === testing
	void crossgrid_internal_tests();
};


#endif
