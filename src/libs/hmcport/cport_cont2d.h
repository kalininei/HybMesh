#ifndef HYBMESH_CPORT_CONT2D_H
#define HYBMESH_CPORT_CONT2D_H
#include "hmcport.h"

//all interface c-functions return HMSUCCESS/HMERROR state values
extern "C"{

int c2_dims(void* obj, int* ret);

int c2_deepcopy(void* obj, void** ret);

int c2_free(void* obj);

int c2_area(void* obj, double* ret);

int c2_length(void* obj, double* ret);

//moves all points using (dp[0], dp[1]) vector
int c2_move(void* obj, double* dp);

//scales all points at pc[0] percent at x, pc[1] percent at y direction
//keeping (p[0], p[1]) point at its place
int c2_scale(void* obj, double* pc, double* p0);

//reflects geometry relative to (v0[0, 1], v1[0, 1]) vector
int c2_reflect(void* obj, double* v0, double* v1);

//rotates object around (p0[0], p0[1]) point at 'a' degrees
int c2_rotate(void* obj, double* p0, double a);

//builds contour using consequtive list of points
//and consequitive list of boundary types.
int c2_frompoints(int npts, double* pts, int* bnds, int force_closed, void** ret);

//simplifies contour so that all resulting sections were at
//no more than *angle* degree. If keep_btypes -> do not merge edges
//with different boundary types.
int c2_simplify_self(void* obj, double angle, int keep_btypes);

//breaks edges graph on independent subgraphs
int c2_quick_separate(void* obj, int* nret, void*** ret);


//unite contours merging equal vertices and edges
int c2_unite(int nobjs, void** objs, void** ret);

//decomposes edges so that resulting contours were either simple open
//or simple closed contours. Collects as many closed contours as possible (even by doubling edges).
int c2_decompose(void* obj, int* nret, void*** ret);


//builds spline which passes through given points and has 'ned' edges
int c2_spline(int np, double* pts, int ned, int* bnds, void** ret);

//clip domain operation: 'union', 'intersection', 'xor', 'difference'
int c2_clip_domain(void* obj1, void* obj2, const char* op, int simplify, void** ret);

//makes contour partition using different algorithms: 'const', 'ref_points', 'ref_weights', 'ref_lengths'
// algo: 0 - const step; 1 - refference point step; 2,3 - ref weights/lengths steps
// n_steps: number of reference points in case of algo=1
// steps: if algo = 'const'      => [const_step],
//        if algo = 'ref_points' => [step0, x0, y0, step1, x1, y1, ...]
//        if algo = 'ref_weights, ref_lengths' => [step0, s0, step1, s1, ...]
// a0: insignificant angle [180-a0, 180 + a0]
// keepbnd: =true if all boundary type changing nodes should be preserved
// n_outbnd - number of output contour edges or -1
int c2_partition(void* obj, const char* algo, int nstep, double* step, double a0, int keepbnd, int nedges,
                 int ncrosses, void** crosses, int nkeeppts, double* keeppts,
                 double* start, double* end,
                 void** ret);

//makes contour partition taking into account partition of given constrained contours and points
int c2_matched_partition(void* obj, int nconts, void** conts, int npts, double* pts, double step,
                         double infdist, double power, double a0, void** ret);

int c2_segment_partition(double start, double end, double hstart, double hend,
	int ninternal, double* hinternal, int* nret, double** ret);

//breaks contour into subcontours by given points
int c2_extract_subcontours(void* obj, int nplist, double* plist, void** ret);

//connects subcontours into single contour using shift and stretch operations
//close_method = 'no', 'yes', 'force'
int c2_connect_subcontours(int nobjs, void** objs, int nfx, int* fx, int shift,
		const char* close_method, void** ret);

//set boundary types to contour
int c2_assign_boundary_types(void* obj, int* bnd, int** revdif);

//-> (0) 'open', (1) 'closed', (2) 'mdomain', (3) 'compound'
int c2_contour_type(void* obj, int* ret);

//assembles a contour by shallow copying primitives of given set
int c2_concatenate(int nobjs, void** objs, void** ret);

//computes closest points using proj = 'vertex'/'edge' algorithms
int c2_closest_points(void* obj, int npts, double* pts, const char* proj, double* ret);

//returns point at given index
int c2_point_at(void* obj, int index, double* ret);

//end points for open contour or throws
int c2_end_points(void* obj, double* start, double* end);

//connectivity tables
int c2_tab_edgevert(void* obj, int* ret);
int c2_tab_vertices(void* obj, double* ret);
int c2_tab_btypes(void* obj, int* ret);

//export operations
int c2_to_hm(void* doc, void* node, void* obj, const char* name, const char* fmt);

}
#endif
