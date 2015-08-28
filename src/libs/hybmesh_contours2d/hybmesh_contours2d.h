#ifndef HYBMESH_CONTOURS2D_H
#define HYBMESH_CONTOURS2D_H

//========== C interface for library
extern "C"{
//Elementary function for linkage testing. Returns 2*a.
int hybmesh_contours2d_ping(int a);

//Contour operations
void  add_point_to_contour(void* cont, double x, double y);
void delete_contour(void*);

//ClosedContour operations
void* create_closed_contour();

//ContourTree Operations
void* create_contour_tree();
void  delete_contour_tree(void*);
//add contour
//returns 0: failed, 1: success
int add_contour_to_tree(void* tree, void* closed_contour);

}




#endif
