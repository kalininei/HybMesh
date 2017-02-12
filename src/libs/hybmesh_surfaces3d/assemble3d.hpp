#ifndef HYBMESH_SURFACE_ASSEMBLE3D_HPP
#define HYBMESH_SURFACE_ASSEMBLE3D_HPP
#include "surface.hpp"

namespace HM3D{ namespace Surface{ namespace Assembler{

//assemble surface from grid.
//makes a shallow copy. do not revert.
FaceData GridSurface(const HM3D::GridData& g);
std::map<int, FaceData> GridSurfaceBType(const HM3D::GridData& g);

//subsurface
FaceData SubSurface(const FaceData& s, const Vertex* v);

//extract smooth surface sections with normal deviations less than angle(deg)
vector<FaceData> ExtractSmooth(const FaceData& s, double angle);

}}

namespace Contour{ namespace Assembler{

//Extracts edges boundary
//!!! all boundary edges should have correct direction
EdgeData ExtractBoundary(const FaceData& a, Vertex v);

//return[0] - positive closed boundaries, return[1] - negative
//!!! faces directions should match
std::array<vector<EdgeData>, 2> ExtractAllBoundaries(const FaceData& a, Vect3 right_normal);

//connect edges at data starting from vertex closest to v till close or end reached.
//chooses direction according to edge which includes v as start vertex
//!!! If edge vertex has more then 2 connections or 
//    multiple edges include v as start point result is undefined (false assert in debug mode)
//!!! Doesn't go backward
EdgeData Connect(const EdgeData& data, Vertex v);



}}
}
#endif
