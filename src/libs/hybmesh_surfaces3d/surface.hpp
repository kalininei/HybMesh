#ifndef SURFACE_GRID3D_HPP
#define SURFACE_GRID3D_HPP
#include "primitives3d.hpp"

namespace HM3D{ namespace Surface{

//retuns false if there are edges with != 2 face connections
bool IsClosed(const FaceData& srf);

//signed volume starting from first three points of first face
//!!! fc should be closed
//!!! all faces should have same direction relative to volume
double Volume(const FaceData& fc);

//assemble surface from grid.
//makes a shallow copy. do not revert.
FaceData GridSurface(const HM3D::GridData& g);
std::map<int, FaceData> GridSurfaceBType(const HM3D::GridData& g);

//subsurface
FaceData SubSurface(const FaceData& s, const Vertex* v);

//extract smooth surface sections with normal deviations less than angle(deg)
vector<FaceData> ExtractSmooth(const FaceData& s, double angle);

//topologically unique rearrange with respect to given edge
//edges and faces data would be reordered, face vector will be permuted.
void FaceRearrange(FaceData& s, const Edge* ed);

//are 'a' and 'b' topologicaly equal
bool MatchTopology(const FaceData& a, const FaceData& b);

//Extracts edges boundary
//!!! all boundary edges should have correct direction
EdgeData ExtractBoundary(const FaceData& a, Vertex v);

//return[0] - positive closed boundaries, return[1] - negative
//!!! faces directions should match
std::array<vector<EdgeData>, 2> ExtractAllBoundaries(const FaceData& a, Vect3 right_normal);

//assign boundary by central vertex and old boundary type.
void SetBoundaryTypes(const FaceData& fvec, std::function<int(Vertex, int)> bfun);


}}

#endif

