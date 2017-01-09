#ifndef HMCONT2D_CONTOUR_HPP
#define HMCONT2D_CONTOUR_HPP
#include "primitives2d.hpp"

namespace HM2D{ namespace Contour{

bool IsContour(const EdgeData&);
bool IsClosed(const EdgeData&);
bool IsOpen(const EdgeData&);

shared_ptr<Vertex> First(const EdgeData&);
shared_ptr<Vertex> Last(const EdgeData&);

//first == last point for closed paths
VertexData OrderedPoints(const EdgeData&);
//like ordered points but removes zero length sections (first == last for closed)
VertexData UniquePoints(const EdgeData&);
//list of all corner points in correct order without doubling
//the last point for closed contours
//for open contours includes first and last points by default
VertexData CornerPoints(const EdgeData&);
//same as corner_points but doubles end points for closed contours
VertexData CornerPoints1(const EdgeData&);

//signed area of closed contour
double Area(const EdgeData&);

//reverses edges vector.
void Reverse(EdgeData&);

//shallow copy 'from' edges to beginning or end of 'to'
//does connection only if 'to' and 'from' have same end point
void Connect(EdgeData& to, const EdgeData& from);

//creates new edge connecting last point and ed
void AddLastPoint(EdgeData& to, shared_ptr<Vertex> ed);

//Return arbitrary inner point for a closed contour
//Direction is not considered
Point InnerPoint(const EdgeData&);

//->(INSIDE, BOUND, OUTSIDE). Direction is ignored.
int WhereIs(const EdgeData&, const Point& p);

// =============== local coordiantes
//gives
//<0>  contour length coordinate, 
//<1>  contour weight coordinate,
//<2>  index of edge
//<3>  edge weight coordinate
//<4>  distance to contour
//of the contour point closest to given one.
std::tuple<double, double, int, double, double>
CoordAt(const EdgeData&, const Point& p);

struct PInfoR{
	shared_ptr<Vertex> p, pprev, pnext;
	shared_ptr<Edge>  eprev, enext;
	int index;
};
//returns vector of length (size()+1)
//detailed information about each node connection.
vector<PInfoR> OrderedInfo(const EdgeData& ed);
PInfoR PInfo(const EdgeData&, const Point*);


Point WeightPoint(const EdgeData&, double);
vector<Point> WeightPoints(const EdgeData&, vector<double>);
vector<Point> WeightPointsByLen(const EdgeData&, vector<double>);

//return weights of points as given by c.ordered_points()
//first is 0, last is always 1.
vector<double> EWeights(const EdgeData&);

//returns true if direction of eind-th edge coinsides with contour direction.
bool CorrectlyDirectedEdge(const EdgeData& dt, int eind);

//Find a point on a edges set, closest to p and place it there by splitting edge.
//Target edge will be removed from input vector, two new edges pasted.
//If p already exists do nothing.
//Returns:
//    <0> if point was placed
//    <1> shared point of point in dt equal to p
std::tuple< bool, shared_ptr<Vertex> >
GuaranteePoint(EdgeData& dt, const Point& p);

}};

#endif
