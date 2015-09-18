#ifndef HYBMESH_HMFEM_FEMGRID43
#define HYBMESH_HMFEM_FEMGRID43
#include "grid.h"

namespace HMFem{ namespace Impl{

//grid shares points with another gridgeom object
//but has connectivity with only 3/4 nodes elements
class Grid43: public GridGeom{
public:
	Grid43(GridGeom* orig);
};


}}
#endif
