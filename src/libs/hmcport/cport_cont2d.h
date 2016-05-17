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
		int* n_outbnd, int** outbnd);


}
#endif
