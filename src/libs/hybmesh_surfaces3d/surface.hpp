#ifndef SURFACE_GRID3D_HPP
#define SURFACE_GRID3D_HPP
#include "primitives3d.hpp"

namespace HM3D{ namespace Surface{

//retuns false if there are edges with != 2 face connections
bool IsClosed(const FaceData& srf);

//signed volume starting from first three points of first face
//!!! fc should be closed
//!!! all faces should have same direction relative to volume
//    use treverter procedure to guarantee this
double Volume(const FaceData& fc);

//topologically unique rearrange with respect to given edge
//edges and faces data would be reordered, face vector will be permuted.
void FaceRearrange(FaceData& s, const Edge* ed);

//are 'a' and 'b' topologicaly equal
bool MatchTopology(const FaceData& a, const FaceData& b);

//assign boundary by central vertex and old boundary type.
void SetBoundaryTypes(const FaceData& fvec, std::function<int(Vertex, int)> bfun);


}}

#endif

