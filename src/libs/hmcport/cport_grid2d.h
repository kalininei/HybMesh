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

//data_periodic is a 
int export_msh_grid(const Grid* grid, const char* fname,
		const Grid2DBoundaryStruct* bstr,
		const BoundaryNamesStruct* bnames,
		int n_periodic,
		int* data_periodic);

}
#endif
