#ifndef HMCONT2D_CONTOUR_HPP
#define HMCONT2D_CONTOUR_HPP
#include "primitives2d.hpp"

namespace HM2D{ namespace Contour{

bool IsContour(const EdgeData&);
bool IsClosed(const EdgeData&);
bool IsOpen(const EdgeData&);

shared_ptr<Vertex> First(const EdgeData&);
shared_ptr<Vertex> Last(const EdgeData&);

//==  Sequences of points for contours
//== first=last point for closed paths
VertexData OrderedPoints(const EdgeData&);
//like ordered points but removes zero length sections
VertexData UniquePoints(const EdgeData&);
//list of all corner points in correct order
//for open contours includes first and last points by default
VertexData CornerPoints(const EdgeData&);
//corner points + points with different left/right boundary types
VertexData SignificantPoints(const EdgeData&);

//==  Sequences without doubling last point for closed contours
VertexData OrderedPoints1(const EdgeData&);
VertexData UniquePoints1(const EdgeData&);
VertexData CornerPoints1(const EdgeData&);
VertexData SignificantPoints1(const EdgeData&);

//returns true if direction of eind-th edge coinsides with contour direction.
bool CorrectlyDirectedEdge(const EdgeData& dt, int eind);

//signed area of closed contour
double Area(const EdgeData&);

//Return arbitrary inner point for a closed contour
//Direction is not considered
Point InnerPoint(const EdgeData&);

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

//lengths
double Length(const EdgeData&);
vector<double> ELengths(const EdgeData&);

//weights
Point WeightPoint(const EdgeData&, double);
vector<Point> WeightPoints(const EdgeData&, vector<double>);
vector<Point> WeightPointsByLen(const EdgeData&, vector<double>);
//get boundary types from contour coordinate
vector<int> BTypesFromWeights(const EdgeData&, const vector<double>&);

//return weights of points as given by c.ordered_points()
//first is 0, last is always 1.
vector<double> EWeights(const EdgeData&);

}};

#endif
