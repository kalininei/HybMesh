#ifndef HMCONT2D_ALGOS_HPP
#define HMCONT2D_ALGOS_HPP

#include "primitives2d.hpp"
#include "contour_tree.hpp"

namespace HM2D{

//algorithms for consequitive contours and
//contour trees
namespace Contour{ namespace Algos{

// ========= Simplifications
//remove points which lie on the same edge
//removes zero length edges
EdgeData Simplified(const EdgeData& cont);
Contour::Tree Simplified(const Tree& t1);

//removes points by indicies
//if inner point is removed then a new edge will be created.
//It will keep all features of the edge previous to deleted point
void RemovePoints(EdgeData& data, vector<int> ipnt);

//Find a point on an edges set, closest to p, and place it there by splitting edge.
//Target edge will be removed from input vector, two new edges pasted.
//If p already exists do nothing.
//Returns:
//    <0> if point was placed
//    <1> shared pointer to point in dt equal to p
std::tuple< bool, shared_ptr<Vertex> >
GuaranteePoint(EdgeData& dt, const Point& p);

//splits edge by addition of sequence of inner points.
//Additions will be done using the direction of the whole contour (not the directin of edge).
//iedge will be conserved and shrinked by the pts[0] point.
//pts.size() new edges will be added
void SplitEdge(EdgeData& cont, int iedge, const vector<Point>& pts);

//reverses edges vector.
void Reverse(EdgeData&);

//shallow copy 'from' edges to beginning or end of 'to'
//does connection only if 'to' and 'from' have same end point
void Connect(EdgeData& to, const EdgeData& from);

//creates new edge connecting last point and ed
void AddLastPoint(EdgeData& to, shared_ptr<Vertex> ed);

}}

// ============================ algorithms for shattered edges collections
namespace ECol{ namespace Algos{
//Does simplifications of simple contours assembled from ecol shattered edges.
//makes deep copies of input edges (not points).
//if angle<0 returns everything back making copies of edges but not points.
EdgeData Simplified(const EdgeData& ecol, double degree_angle, bool no_break_at_bt_change=false,
		const VertexData& keep={});

//deep copies everything
EdgeData NoCrosses(const EdgeData& ecol);

void MergePoints(EdgeData& ecol);

//set 'to' boundary types from closest 'from' edges.
void AssignBTypes(const EdgeData& from, EdgeData& to);
}}


}
#endif
