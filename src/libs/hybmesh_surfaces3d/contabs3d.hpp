#ifndef HYBMESH_CONTABS3D_HPP
#define HYBMESH_CONTABS3D_HPP

#include "primitives3d.hpp"

namespace HM3D{ namespace Connectivity {

//Vertex->Edges
struct VertexEdgeR{
	shared_ptr<Vertex> v;
	vector<int> eind;
	size_t size() const { return eind.size(); }
};
vector<VertexEdgeR> VertexEdge(const ShpVector<Edge>& data);

//Edge->Face
struct EdgeFaceR{
	shared_ptr<Edge> e;
	vector<int> find;  //faces indicies
	size_t size() const { return find.size(); }
};
vector<EdgeFaceR> EdgeFace(const FaceData& data);

//Edge->Face extended
//   find - face index,
//   locind - local edge index within the face,
//   posdir - is edge directed according to face
struct EdgeFaceExtendedR: public EdgeFaceR{
	vector<int> locind;
	vector<bool> posdir;
};
vector<EdgeFaceExtendedR> EdgeFaceExtended(const FaceData& data);

//Face->Face
vector<vector<int>> FaceFace(const FaceData& data);
vector<vector<int>> FaceFace(const vector<EdgeFaceR>& edge_face, int nfaces);

//Cell->Cell
vector<vector<int>> CellCell(const CellData& data);


}}


#endif
