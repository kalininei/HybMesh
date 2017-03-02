#ifndef HYBMESH_HMCPORT_GRID2D_H
#define HYBMESH_HMCPORT_GRID2D_H

#include "hmcport.h"

extern "C"{

int g2_dims(void* obj, int* ret);
int g2_area(void* obj, double* ret);
int g2_bnd_dims(void* obj, int* ret);
int g2_bnd_length(void* obj, double** ret);
int g2_skewness(void* obj, double threshold, double* maxskew, int* maxskewindex,
		int* badnum, int** badindex, double** badvals);
int g2_deepcopy(void* obj, void** ret);
int g2_free(void* obj);
int g2_concatenate(int nobjs, void** objs, void** ret);
int g2_move(void* obj, double* dx);
int g2_scale(void* obj, double* pc, double* p0);
int g2_reflect(void* obj, double* v0, double* v1);
int g2_rotate(void* obj, double* p0, double a);


int g2_from_points_edges(int npoints, double* points, int neds, int* eds, void** ret);
int g2_from_points_cells(int npoints, double* points, int ncells, int* cellsizes, int* cellvert,
		int nbedges, int* bedges, void** ret);
int g2_rect_grid(int nx, double* xdata, int ny, double* ydata, int* bnds, void** ret);
int g2_circ_grid(double* p0, int nr, double* rdata, int na, double* adata, int istrian, int bnd, void** ret);
int g2_ring_grid(double* p0, int nr, double* rdata, int na, double* adata, int* bnds, void** ret);
int g2_tri_grid(double* verts, int nedge, int* bnds, void** ret);
//areatype:
//  'rect' -- rectangle given by two points
//  'hex' -- hexagon given by center point and radius
int g2_hex_grid(const char* areatype, double* area, double crad, int strict, void** ret);
int g2_extract_contour(void* obj, void** ret);
int g2_assign_boundary_types(void* obj, int* bnd, int** revdif);

//triangulate given contour using contour vertices as basic points
// domain - domain to triangulate
// constr - constraint or NULL
// nemb - number of embedded points
// emb - set of data in [x1, y2, size2, x2, y2, size2, ...] representing embedded points
// filler = '3', '4', 'pebi'
//returns pointer to HM2D::GridData or NULL if failed
int g2_unstructured_fill(void* domain, void* constraint,
		int nembpts, double* embpts, const char* filler, void** ret);

//builds rectangular grid on basis of four open HMCont2D::Contour void* objects
//algo is:
//  0) linear
//  1) inverse_laplace
//  2) direct_laplace
//  3) orthogonal
//  4) linear_tfi
//  5) hermite_tfi
//herw: dx weights for left, bot, right, top.
//   dx=1 => normal connection, dx=0 => linear connection
//if input are not positioned incorrectly their vertices coordinates will be changed
int g2_custom_rect_grid(const char* algo, void* left, void* bottom, void* right, void* top,
		double* herw, int rinvalid, void** ret, hmcport_callback cb);

int g2_circ4grid(const char* algo, double* p0, double rad, double step, double sqrside, double rcoef, void** ret);

int g2_stripe_grid(void* obj, int npartition, double* partition, const char* tipalgo, int* bnd, void** ret,
		hmcport_callback cb);

//base_grid - Grid* object 
//target_contour - EdgeData* object
//pbase, ptarget[2*Npnt] - mapped points given as [x0, y0, x1, y1, ....]
//snap_method:
//   * - no
//   * - add_vertices
//   * - shift_vertices
//algo:
//   1 - direct_laplace
//   2 - inverse_laplace
//reversed:
//   0 - default
//   1 - target contour is reversed
//if return_invalid = true then grid will be returned even if it is not valid,
int g2_map_grid(void* base_obj, void* target_obj,
		int npoints, double* base_points, double* target_points,
		const char* snap, int bt_from_contour, const char* algo,
		int is_reversed, int rinvalid, void** ret, hmcport_callback cb);

int g2_point_at(void* obj, int index, double* ret);
int g2_closest_points(void* obj, int npts, double* pts, const char* proj, double* ret);
int g2_simplify_bnd(void* obj, double angle, void** ret);


//tabs
int g2_tab_btypes(void* obj, int* ret);
int g2_tab_vertices(void* obj, double* ret);
int g2_tab_edgevert(void* obj, int* ret);
int g2_tab_cellsizes(void* obj, int* ret);
int g2_tab_cellvert(void* obj, int* nret, int** ret2);
int g2_tab_celledge(void* obj, int* nret, int** ret2);
int g2_tab_centers(void* obj, double* ret);
int g2_tab_bedges(void* obj, int* nret, int** ret2);
int g2_tab_edgecell(void* obj, int* ret);
int g2_tab_bndbt(void* obj, int* nret, int** ret2);

//angle in (0, 180). If 0 - removes all concave cell segments. If 180 - only degenerate
//return 0 on error
int g2_convex_cells(void* obj, double angle, void** ret);
int g2_exclude_cont(void* obj, void* cont, int isinner, void** ret, hmcport_callback cb);
int g2_unite_grids(void* obj1, void* obj2, double buf, int fixbnd,
		int emptyholes, double angle0, const char* filler,
		void** ret, hmcport_callback cb);
int g2_to_msh(void* obj, const char* fname, BoundaryNamesStruct btypes, int n_per_data, int* per_data);
int g2_to_tecplot(void* obj, const char* fname, BoundaryNamesStruct btypes);
int g2_to_hm(void* doc, void* node, void* obj, const char* name, const char* fmt, int naf, const char** af);

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
int g2_boundary_layer(int nopt, BoundaryLayerGridOption* opt, void** ret, hmcport_callback cb);


/*
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

//builds rectangular grid on basis of four open HMCont2D::Contour void* objects
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
//returns pointer to HM2D::GridData or NULL if failed
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
//returns pointer to HM2D::GridData or NULL if failed
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
*/

}

#endif
