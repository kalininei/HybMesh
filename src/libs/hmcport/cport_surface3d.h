#ifndef HYBMESH_HMCPORT_SURFACE3D_H
#define HYBMESH_HMCPORT_SURFACE3D_H

extern "C"{

typedef void CPortSurface3D;  //SSurface*

void free_srf3(CPortSurface3D*);

void surf3_dims(const CPortSurface3D*, int* dims);

//returns 0 on faile and 1 on success
//fills ret array 
int surf3_btypes(const CPortSurface3D*, int* ret);

//returns 0 on fail or SSurface pointer on success
void* extract_grid3_surface(const void* grid3);

//return 0 on fail ant 1 on success
//fills (*ret) with (*num) pointers to SSurface
int extract_subsurfaces(const CPortSurface3D*, int* num, CPortSurface3D*** ret);

//returns not null on succes
void* s3writer_create(const char* sname, CPortSurface3D* surf, void* awriter,
		void* subnode, const char* fmt);
int s3writer_add_defined_field(void* swriter, const char* field);
void s3writer_free(void* swriter);

void* s3reader_create(void* awriter, void* subnode, char* outname);
CPortSurface3D* s3reader_getresult(void* rd);
void s3reader_free(void* sreader);

}
#endif

