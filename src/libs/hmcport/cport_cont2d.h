#ifndef HYBMESH_CPORT_CONT2D_H
#define HYBMESH_CPORT_CONT2D_H
#include "hmcport.h"

extern "C"{


// ==================== Contour Partition algorithm
// cont: ECOllection pointer
// btypes: size = cont.size(). Boundary feature for each contour edge
// algo: 0 - const step; 1 - refference point step
// n_steps: number of reference points in case of algo=1
// steps: if algo = 0 => [const_step], if algo = 1 => [step0, x0, y0, step1, x1, y1, ...]
// a0: insignificant angle [180-a0, 180 + a0]
// keepbnd: =true if all boundary type changing nodes should be preserved
// n_outbnd - number of output contour edges
// outbnd - boundary feature for each output contour edge
// nedges - required number of edges or -1 if it should be computed from steps
// returns new ECollection pointer or NULL if failed
void* contour_partition(void* cont, int* btypes, int algo,
		int n_steps, double* steps, double a0, int keepbnd, int nedges,
		int n_crosses, void** crosses,
		int* n_outbnd, int** outbnd);

//build a spline using basis points
void* spline(int npnt, double* pnts, int nbtypes, int* btypes, int nedges,
		int* n_outbnd, int** outbnd);

// builds a partition of contour 'cont' using constant 'step'
// and existing parititions of conds[i] contours.
void* matched_partition(void* cont, int ncond, void** conds, int npts, double* pcond, double step,
		double influence_dist, double pw, double a0);

//returns 0 if failed
int segment_part(double start, double end, double h0, double h1,
		int n_internals, double* h_internals, int* nout, double** hout);


//projection methods: line, vertex, corner
void* extract_contour(void* source, double* pnts, const char* method);

//return 0 on fail
int unite_contours(int ncont, void** conts, void** retcont, int* Nlinks, int** links);
int simplify_contour(void* cont, double degree_angle, int* btypes, void** ret_cont, int* Nretb, int** retb);
int separate_contour(void* cont, int* btypes, int* Nretc, void*** retc, int* Nretb, int** retb);
int quick_separate_contour(void* cont, int* btypes, int* Nretc, void*** retc, int* Nretb, int** retb);

void* cwriter_create(const char* cname, void* cont, void* awriter, void* subnode, const char* fmt);
void cwriter_free(void* cwriter);
int cwriter_add_edge_field(void* cwriter, const char* fieldname, void* field, int fsize, const char* type);
void* creader_create(void* awriter, void* subnode, char* outname);
void* creader_getresult(void* rd);
void* creader_read_edge_field(void* rd, const char* fieldname, const char* type);
void creader_free(void* creader);

}
#endif
