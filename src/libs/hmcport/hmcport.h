#ifndef HYBMESH_HMCPORT_H
#define HYBMESH_HMCPORT_H

#include "crossgrid.h"

extern "C"{


// === grid management
Grid* grid_construct(int Npts, int Ncells, double* pts, int* cells);
//grid structures count
int grid_npoints(Grid* g);
int grid_ncells(Grid* g);
int grid_cellsdim(Grid* g);
//area
double grid_area(Grid* g);
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
//the same without preallocation
void grid_get_points_cells2(Grid* g, double** pts, int** cells);

//get edge->cell connectivity table: length of resulting ed_cell is 2*Neds;
//ed_pt is edge->points connectivity which is used for sorting of the resulting edges.
//if ed_pt=NULL -> resulting edges will be sorted arbitraty
//returns 1 if ok or 0 if error input (edge->points table is invalid).
int grid_get_edge_cells(Grid* g, int* Neds, int** ed_cell, int* ed_pt);
void grid_free_edge_cells(int** ed_cell);

//save grid to vtk
void grid_save_vtk(Grid* g, const char* fn);
//free grid
void grid_free(Grid* g);

// === contours management
//Creates a contour: from edges->points connectivity
//Input data should be assembled as a list of closed polygons
//with inner first level contour, outer second level etc
//inner ones have an anti-clockwise direction, outer - clockwise.
//Returns NULL if input data are not assembled correctly
Cont* contour_construct(int Npts, int Ned, double* pts, int* edges);
//Frees the contour
void cont_free(Cont* c);
//save contour to vtk
void contour_save_vtk(Cont* c, const char* fn);

//get contours points and edges->points connectivity
void contour_get_info(Cont* c, int* Npnt, int* Neds, 
		double** pts,
		int** eds);
//free edges data
void contour_free_info(double** pts, int** eds);
//area
double contour_area(Cont* c);

//Fill target boundary values from source contour values if
//target edge lies on source edge. Else assigns def value.
void add_contour_bc(Cont* src, Cont* tar, int* vsrc, int* vtar, int def);

// === hybmesh_contours2d management
//builds HMBlay::Container<HMBlay::ECollection>.
//pts[2*Npts], edges[2*Nedges]
void* create_ecollection_container(int Npts, double* pts, int Nedgs, int* edges);
void free_ecollection_container(void* ecol);
//calculates area with respect to multiplicity
double ecollection_area(void* ecol);

// === set default callbacks
//crossgrid callback pointer is the function of the following arguments:
//	(global procedure name, local procedure name, global percentage, local percentage)
//percentages are doubles in [0, 1] range. If it is less then, then percentages of procedures is not tracking.
//normally returns CrossGridCallback::OK. Should return CrossGridCallback::CANCEL for cancellation requiry.
typedef int (*crossgrid_callback)(const char*, const char*, double, double);
void crossgrid_cout_callback();  //callback to std::cout
void crossgrid_silent_callback(); //no callback at all
void crossgrid_set_callback(crossgrid_callback fun); //user defined callback

// === crossgrid procedure
// buffer_size -- distance from gsecondary contour which defines buffer zone.
//                Could be zero
// density -- relative density of triangles in a buffer zone in [0,1]
// preserve_bp -- equals to 1 if all boundary nodes should present in resulting grid.
//                If 0 some non-corner nodes could be eliminated or shifted
// empty_holes -- equals to 1 if all holes in secondary grid should present in
//                resulting grid. Otherwise (=0) these holes will be filled with gbase mesh.
//with callback defined globally
Grid* cross_grids(Grid* gbase, Grid* gsecondary, double buffer_size, 
		int preserve_bp, int empty_holes);
//with specified callback function
Grid* cross_grids_wcb(Grid* gbase, Grid* gsecondary, double buffer_size, 
		int preserve_bp, int empty_holes, crossgrid_callback cb_fun);

// === NewGrid = Grid exclude Contour Area
//with callback defined globally
Grid* grid_exclude_cont(Grid* grd, Cont* cont, int is_inner);
//with specified callback function
Grid* grid_exclude_cont_wcb(Grid* grd, Cont* cont, int is_inner,
		crossgrid_callback cb_fun);


// === Boundary layer grid constructing
// Build a boundary layer grid around a contour tree
// Returns 0 if fails or resulting Grid object.
struct BoundaryLayerGridOption{
	void* cont;  // HMCont2::ECollection
	int Npart;  //length of part
	double* part; //partition array starts with 0
	const char* tp;  //INSIDE, OUTSIDE, LEFT, RIGHT, AROUND
	const char* mesh_cont; //'NO', 'KEEP_ORIGIN', 'KEEP_SHAPE', 'IGNORE_ALL'
	double mesh_cont_step; //contour step if mesh_cont is not 'NO'
	double start[2]; //start point
	double end[2];   //end point
	int force_conformal;   //force conformal
	double angle_range[4]; //maximum acute/right/straight/reentrant angles
};

Grid* boundary_layer_grid_wcb(int N, BoundaryLayerGridOption* opt, 
		crossgrid_callback cb_fun);

}; //extern C


#endif
