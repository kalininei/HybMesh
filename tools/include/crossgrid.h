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
	int grid_npoints(Grid* g);
	int grid_ncells(Grid* g);
	int grid_cellsdim(Grid* g);
	void grid_get_points_cells(Grid* g, double* pts, int* cells);
	void grid_save_vtk(Grid* g, const char* fn);
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
