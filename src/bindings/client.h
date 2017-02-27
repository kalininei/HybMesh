#ifndef CORE_HYBMESH_CONNECTION_H
#define CORE_HYBMESH_CONNECTION_H

#include "stdint.h"

#ifdef __cplusplus
extern "C"{
#endif

typedef struct _hybmesh_client2server{
	int isopen;
	int sig_read;
	int sig_write;
	int data_read;
	int data_write;
	intptr_t childid;
} HybmeshClientToServer;


/* return 0 on errors */
int HybmeshClientToServer_new(const char* path, int* id);
int HybmeshClientToServer_get_signal(int id, char* sig);
int HybmeshClientToServer_get_data(int id, int* sz, char** data);
int HybmeshClientToServer_get_data1(int id, int* sz);
int HybmeshClientToServer_get_data2(int id, int sz, char* data);
int HybmeshClientToServer_send_signal(int id, char sig);
int HybmeshClientToServer_send_data(int id, int sz, char* data);
int HybmeshClientToServer_delete(int id);

#ifdef __cplusplus
}
#endif

#endif
