#ifndef HYBMESH_HMCPORT_H
#define HYBMESH_HMCPORT_H

#include <stddef.h>

#define HMSUCCESS 1
#define HMERROR 0

//all interface c-functions return HMSUCCESS/HMERROR state values
extern "C"{

//callback arguments: procedure name,
//                    subprocedure name,
//                    [0, 1] procedure done,
//                    [0, 1] subprocedure done
//         returns:   HMCallback::OK (=0) to proceed with operation
//                    HMCallback::Cancel (=1) to cancel operation
typedef int (*hmcport_callback)(const char*, const char*, double, double);

//calculates ascii file hash
//all floating points in file will be rounded to 0.00 before calculating hash
//so file with content "0.711 23.00000001" would have same hash as "0.71 23"
int get_ascii_file_hash(const char* fn, size_t* ret);

//frees array allocated on c-side by new[] commands
int free_int_array(int* a);
int free_double_array(double* a);
int free_char_array(char* a);
int free_voidp_array(void** a);

//returns last error message passed using add_error_message
//by any of hmcport functions.
int get_last_error_message(char** msg);
void add_error_message(const char*);

// structure representing boundary-index->boundary-name dictionary
struct BoundaryNamesStruct{
	int n;
	int* index;
	const char** name;
};
const char* get_boundary_name(const BoundaryNamesStruct*, int index);
}; //extern C


#endif
