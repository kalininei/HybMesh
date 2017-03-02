#ifndef HYBMESH_HMCPORT_GRID3D_H
#define HYBMESH_HMCPORT_GRID3D_H

#include "hmcport.h"

extern "C"{

int g3_move(void* obj, double* dx);
int g3_scale(void* obj, double* pc, double* p0);
int g3_point_at(void* obj, int index, double* ret);

int g3_deepcopy(void* obj, void** ret);
int g3_concatenate(int nobjs, void** objs, void** ret);

//====== destructor
int g3_free(void* obj);

//====== information
//n_vert, n_edges, n_faces, n_cells
int g3_dims(void* obj, int*);

//n_vert, n_edges, n_faces
int g3_bnd_dims(void* obj, int* dims);

//fills ret array with 'number of faces' values
int g3_tab_btypes(void* obj, int* ret);

int g3_tab_vertices(void* obj, double* ret);
int g3_tab_edgevert(void* obj, int* ret);
int g3_tab_facedim(void* obj, int* ret);
int g3_tab_faceedge(void* obj, int* nret, int** ret2);
int g3_tab_facevert(void* obj, int* nret, int** ret2);
int g3_tab_facecell(void* obj, int* ret);
int g3_tab_cellfdim(void* obj, int* ret);
int g3_tab_cellvdim(void* obj, int* ret);
int g3_tab_cellface(void* obj, int* nret, int** ret2);
int g3_tab_cellvert(void* obj, int* nret, int** ret2);
int g3_tab_bnd(void* obj, int* nret, int** ret2);
int g3_tab_bndbt(void* obj, int* nret, int** ret2);

//boundary area
int g3_bnd_area(void* obj, double* ret);

//creates surface out of grid boundary
int g3_extract_surface(void* obj, void** ret);

//volume
int g3_volume(void* obj, double* ret);

//merge coincident primitives
int g3_merge(void* obj1, void* obj2, void** ret, hmcport_callback cb);

int g3_assign_boundary_types(void* obj, int* bnd, int** revdif);

//====== sweep constructors
//construct by sweep in z direction
//btop, bbot = boundary types for each 2d cell of btop and bbot boundaries
//bside value for side boundary (-1 to take values from 2d grid)
int g3_extrude(void* obj, int nz, double* zvals,
		int*  bbot, int* btop,
		int bside, void** ret);

//vec - [x0, y0, x1, y1] array defining vector of rotation
//phi[n_phi] - increasing vector of angular partition (degree)
//b1, b2 - boundary types for surfaces at minimum and maximum phi's
//is_trian (bool) - whether to triangulate center cell
//return NULL if failed
int g3_revolve(void* obj, double* vec, int n_phi, double* phi,
		int is_trian, int b1, int b2, void** ret);


//======= unstructured fill
//0 on fail
int g3_tetrahedral_fill(int nsurf, void** surf,
		int nconstr, void** constr,
		int npts, double* pcoords, double* psizes,
		void** ret, hmcport_callback cb);


//====== exporters
int g3_to_vtk(void* obj, const char* fname, hmcport_callback f2);
int g3_surface_to_vtk(void* obj, const char* fname, hmcport_callback f2);
int g3_to_msh(void* obj, const char* fname, BoundaryNamesStruct bnames,
		int n_periodic, double* data_periodic, hmcport_callback f2);
int g3_to_gmsh(void* obj, const char* fname, BoundaryNamesStruct bnames,
		hmcport_callback f2);
int g3_to_tecplot(void* obj, const char* fname, BoundaryNamesStruct bnames,
		hmcport_callback f2);
int g3_to_hm(void* doc, void* node, void* obj, const char* name, const char* fmt,
		int naf, const char** af,
		hmcport_callback f2);
}

#endif
