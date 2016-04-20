#ifndef HYBMESH_HMCPORT_GRID2D_H
#define HYBMESH_HMCPORT_GRID2D_H
#include "hmcport.h"

extern "C"{

//boundary struct
struct Grid2DBoundaryStruct{
	int n;
	int* edge_start_nodes;
	int* edge_end_nodes;
	int* btypes;
};
Grid2DBoundaryStruct* set_grid2_boundary_types(int n, int* n1, int* n2, int* bt);
void free_grid2_boundary_types(Grid2DBoundaryStruct*);

//data_periodic is a an array of 3*n_periodic ints like
//	[boundary-per, boundary-shadow, is-reversed, ...]
//returns 0 on success
int export_msh_grid(const Grid* grid, const char* fname,
		const Grid2DBoundaryStruct* bstr,
		const BoundaryNamesStruct* bnames,
		int n_periodic,
		int* data_periodic);

int export_tecplot_grid(const Grid* grid, const char* fname,
		const Grid2DBoundaryStruct* bstr,
		const BoundaryNamesStruct* bnames);

//builds rectangular grid on basis of four open HMCont2D::Contour objects
//algo is:
//  0) linear connection
//  1) laplas connection
//  2) conformal connection
//if input are not positioned correctly their vertices coordinates will be changed
//returns pointer to GridGeom or NULL if failed
void* custom_rectangular_grid(int algo, void* left, void* bot, void* right, void* top);

//builds a grid with quadrangle cells in a circular area. it contains uniform square grid
//in the center and ring-like grid at the outer boundary.
//algo:
//	0: laplas8
//center[2]: x, y of center point
//rad: radius of circle area
//step: step size at the outer boundary
//sqrside: side of the inner square normalized by radius (1.0 is good, >1.4 is error)
//outer_refinement: radius step normalized by arc step at the outer boundary (1 - no refinement, < 1 - with refinement)
//returns: NULL if fails
Grid* circ4grid(int algo, double* center, double rad, double step, double sqrside, double outer_refinement);




}
#endif
