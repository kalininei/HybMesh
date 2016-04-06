#ifndef HYBMESH_HMCPORT_H
#define HYBMESH_HMCPORT_H

#include "crossgrid.h"
#include "hmcallback.hpp"

extern "C"{


// === grid management
Grid* grid_construct(int Npts, int Ncells, double* pts, int* cells);
//grid structures count
int grid_npoints(Grid* g);
int grid_ncells(Grid* g);
int grid_cellsdim(Grid* g);
//area
double grid_area(Grid* g);
//equiangular skewness
//returns 0 if ok and 1 on errors
int report_skewness(void* grid, double threshold, double* max_skew, int* max_skew_cell,
		int* bad_cells_count, int* bad_indicies, double* bad_skew);


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

//get ecollection structure
void ecollection_edges_info(void* ecol, int* npts, int* neds, double** pts, int** eds);
void free_ecollection_edges_info(double* pts, int* eds);

//calculates area with respect to multiplicity
double ecollection_area(void* ecol);

// c1, c2 - ecollection pointers
// oper = 1 (union), 2 (difference), 3 (intersection), 4 (xor)
// simplify 1 (true) 0 (false)
// Returns Null if result is empty or error occurred
void* domain_clip(void* c1, void* c2, int oper, int simplify);

//Fill target boundary values from source contour values if
//target edge lies on source edge.
//def - value which will be assigned if no source boundary was found
//src, tar - ecollection objects for source and target boundaries
//vsrc[src.size()], vtar[tar.size()] - input and output boudnary values
//return 0 if ok and 1 if error
int set_ecollection_bc(void* src, void* tar, int def, int* vsrc, int* vtar);

//same as previous but forces boundary assignment depending on algo feature
// 1 - assign only edges which lie on source
// 2 - assign edges which endpoints lie on source
// 3 - assign all
//vsrc is filled only for those values which were found.
int set_ecollection_bc_force(void* src, void* tar, int* vsrc, int* vtar, int algo);


// === crossgrid procedure
// buffer_size -- distance from gsecondary contour which defines buffer zone.
//                Could be zero
// density -- relative density of triangles in a buffer zone in [0,1]
// preserve_bp -- equals to 1 if all boundary nodes should present in resulting grid.
//                If 0 some non-corner nodes could be eliminated or shifted
// empty_holes -- equals to 1 if all holes in secondary grid should present in
//                resulting grid. Otherwise (=0) these holes will be filled with gbase mesh.
// angle0 (degree) - insignificant deviation from the straight angle
//with callback defined globally
Grid* cross_grids(Grid* gbase, Grid* gsecondary, double buffer_size, 
		int preserve_bp, int empty_holes, double angle0);
//with specified callback function
Grid* cross_grids_wcb(Grid* gbase, Grid* gsecondary, double buffer_size, 
		int preserve_bp, int empty_holes, double angle0, HMCallback::Fun2 cb_fun);

// === NewGrid = Grid exclude Contour Area
//with callback defined globally
Grid* grid_exclude_cont(Grid* grd, Cont* cont, int is_inner);
//with specified callback function
Grid* grid_exclude_cont_wcb(Grid* grd, Cont* cont, int is_inner,
		HMCallback::Fun2 cb_fun);

//merges boundary edges of the same cell with angle less than 'angle' in [0, 180]
//return 0 if ok
int simplify_grid_boundary(Grid* grd, double angle);

// === Boundary layer grid constructing
// Build a boundary layer grid around a contour tree
// Returns 0 if fails or resulting Grid object.
struct BoundaryLayerGridOption{
	void* cont;  // HMCont2::ECollection
	int Npart;  //length of part
	double* part; //partition array starts with 0
	const char* tp;  //INSIDE, OUTSIDE, LEFT, RIGHT, AROUND
	const char* mesh_cont; //'NO', 'KEEP_ORIGIN', 'KEEP_SHAPE', 'IGNORE_ALL', 'INCREMENTAL'
	double mesh_cont_step; //contour step if mesh_cont is not 'NO'
	double start[2]; //start point
	double end[2];   //end point
	int force_conformal;   //force conformal
	double angle_range[4]; //maximum acute/right/straight/reentrant angles (deg)
	double step_start; //step values for incremental stepping
	double step_end;
};

Grid* boundary_layer_grid_wcb(int N, BoundaryLayerGridOption* opt, 
		HMCallback::Fun2 cb_fun);


// === Mapping
//base_grid - Grid* object 
//target_contour - ECollection* object
//pbase, ptarget[2*Npnt] - mapped points given as [x0, y0, x1, y1, ....]
//snap_method:
//   1 - no snapping
//   2 - snap with addition of new nodes
//   3 - snap by shifting points without addition new vertices
Grid* build_grid_mapping(void* base_grid, void* target_contour, int Npnt, double* pbase, double* ptarget, int snap_method);


}; //extern C


#endif
