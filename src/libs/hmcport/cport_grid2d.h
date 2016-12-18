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
//  1) inverse-laplace connection
//  2) direct-laplace connection
//  3) orthogonal connection
//  4) linear-tfi
//  5) hermite-tfi
//tfi_hermite_w: dx weights for left, bot, right, top.
//   dx=1 => normal connection, dx=0 => linear connection
//if input are not positioned correctly their vertices coordinates will be changed
//returns pointer to GridGeom or NULL if failed
void* custom_rectangular_grid(int algo, void* left, void* bot,
		void* right, void* top, double* tfi_hermite_w, int return_invalid, hmcport_callback cb_fun);

//builds a grid with quadrangle cells in a circular area. it contains uniform square grid
//in the center and ring-like grid at the outer boundary.
//algo:
//	0: linear
//	1: laplace
//	2: orthogonal-circ
//	3: orthogonal-rect
//center[2]: x, y of center point
//rad: radius of circle area
//step: step size at the outer boundary
//sqrside: side of the inner square normalized by radius (1.0 is good, >1.4 is error)
//outer_refinement: radius step normalized by arc step at the outer boundary (1 - no refinement, < 1 - with refinement)
//returns: NULL if fails
Grid* circ4grid(int algo, double* center, double rad, double step, double sqrside, double outer_refinement);

//triangulate given contour using contour vertices as basic points
// domain - domain to triangulate
// constr - constraints or void
// nemb - number of embedded points
// emb - set of data in [x1, y2, size2, x2, y2, size2, ...] representing embedded points
//algo = 0 -> triangulation
//algo = 1 -> quadrangulation
//returns pointer to GridGeom or NULL if failed
void* triangulate_domain(void* domain, void* constr, int nemb, double* emb, int algo);

void* pebi_fill(void* domain, void* constr, int nemb, double* emb);

//an - angle (0, 180). If 0 - removes all concave cell segments. If 180 - only degenerate
//return 0 on error
void* convex_cells(void* input_grid, double an);

//tip_algo:
//  0 - no tip grid
//  1 - radial tip grid
void* stripe_grid(void* input_contour, int npart, double* part, int tip_algo, 
		void** bot, void** left, void** top, void** right, hmcport_callback cb);


//area_type:
//  =0 -- rectangle given by two points
//  =1 -- hexagon given by center point and radius
void* regular_hex_grid(double* area, int area_type, double cell_rad, int strict_area);


// =========================== writing procedures
//awriter -- ReaderA instance
//subnode -- Reader* instance which is a child node of ReaderA document
//fmt -- "ascii", "bin", "fbin"
//returns Reader* or 0 if error
void* gwriter_create(const char* gname, void* grid, void* awriter, void* subnode, const char* fmt);
void gwriter_free(void* gwriter);
int gwriter_add_defined_field(void* gwriter, const char* field);
int gwriter_add_edge_field(void* gwriter, const char* fieldname, void* field, int fsize, const char* type);
void* greader_create(void* awriter, void* subnode, char* outname);
void* greader_getresult(void* rd);
void* greader_read_edge_field(void* rd, const char* fieldname, const char* type);
void greader_free(void* greader);

}

#endif
