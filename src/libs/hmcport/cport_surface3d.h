#ifndef HYBMESH_HMCPORT_SURFACE3D_H
#define HYBMESH_HMCPORT_SURFACE3D_H

extern "C"{

typedef void CPortSurface3D;  //SSurface*

//returns 0 on fail or SSurface pointer on success
void* extract_grid3_surface(const void* grid3);

}
#endif

