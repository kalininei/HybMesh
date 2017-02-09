#ifndef HYBMESH_HMCPORT_HMXML_H
#define HYBMESH_HMCPORT_HMXML_H

#include "hmcport.h"

//all interface c-functions return HMSUCCESS/HMERROR state values
extern "C"{

int hmxml_open_doc(const char* fname, void** doc, void** root);

int hmxml_new_doc(void** doc, void** root);

int hmxml_free_node(void* node);

int hmxml_free_doc(void* doc);

int hmxml_query(void* node, const char* query, int* rnum, void*** ret);

int hmxml_write(void* doc, const char* filename);

int hmxml_purged_string(void* doc, char** ret);

int read_contour2(void* doc, void* node, void** obj, char* name);

int read_grid2(void* doc, void* node, void** obj, char* name);

int read_surface3(void* doc, void* node, void** obj, char* name);

int read_grid3(void* doc, void* node, void** obj, char* name, hmcport_callback cb);

}
#endif

