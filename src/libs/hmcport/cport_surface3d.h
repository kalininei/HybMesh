#ifndef HYBMESH_HMCPORT_SURFACE3D_H
#define HYBMESH_HMCPORT_SURFACE3D_H
#include "hmcport.h"

extern "C"{

int s3_tab_btypes(void* obj, int* ret);
int s3_tab_vertices(void* obj, double* ret);
int s3_tab_edgevert(void* obj, int* ret);
int s3_tab_facedim(void* obj, int* ret);
int s3_tab_faceedge(void* obj, int* nret, int** cret);
int s3_tab_facevert(void* obj, int* nret, int** cret);
int s3_tab_centers(void* obj, double* ret);


int s3_dims(void* obj, int* dims);
int s3_area(void* obj, double* ret);
int s3_deepcopy(void* obj, void** ret);

int s3_free(void* obj);
int s3_volume(void* obj, double* ret);
int s3_move(void* obj, double* dx);
int s3_scale(void* obj, double* pc, double* p0);
int s3_quick_separate(void* obj, int* nret, void*** ret);
int s3_concatenate(int nobjs, void** objs, void** ret);
int s3_to_hm(void* doc, void* node, void* obj, const char* name, const char* fmt,
		hmcport_callback cb);
int s3_assign_boundary_types(void* obj, int* bnd, int** revdif);

}
#endif

