#ifndef HYBMESH_CONTOURS2D_H
#define HYBMESH_CONTOURS2D_H

//========== C interface for library
extern "C"{
//Elementary function for linkage testing. Returns 2*a.
int hybmesh_contours2d_ping(int a);

//builds a tree of closed contours from a set of connected points
// Npnt - number of points
// points - [2*Npnt] (x0, y0, x1, y1, ...) array of points coords
// Nedges - total number of edges
// edges - [3*Nedges](i0-start, i0-end, i0-btype, i1-start, ...)
//         points to edges + boundary type of each edge 
// Returns ContourTree object or 0 if fails
void* create_contour_tree(int Npnt, double* points,
		int Nedges, int* edges);


//free allocated memory
void free_contour_tree(void* tree);

}
#endif
