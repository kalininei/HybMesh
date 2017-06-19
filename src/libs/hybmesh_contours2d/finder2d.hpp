#ifndef HYBMESH_FINDER2D_HPP
#define HYBMESH_FINDER2D_HPP
#include "primitives2d.hpp"
#include "contabs2d.hpp"
#include "contour_tree.hpp"

namespace HM2D{ namespace Finder{

//contains procedures: returns data object or nullptr
shared_ptr<Vertex> Contains(const VertexData& data, const Point* pnt);
shared_ptr<Vertex> Contains(const EdgeData& data, const Point* pnt);
shared_ptr<Edge> Contains(const EdgeData& data, const Edge* ed);
shared_ptr<Vertex> Contains(const CellData& data, const Point* pnt);
shared_ptr<Edge> Contains(const CellData& data, const Edge* ed);
shared_ptr<Cell> Contains(const CellData& data, const Cell* c);

//closest edge-> returns 
//<0> edge index, -1 if dt is empty
//<1> distance,
//<2> weight of closest point within edge.
std::tuple<int, double, double>
ClosestEdge(const EdgeData& dt, const Point& p);

//closest point lying on edge
Point ClosestEPoint(const EdgeData& dt, const Point& p);

//pointer to closest point:
//<0> - point index
//<1> - squared distance to point
std::tuple<int, double> ClosestPoint(const VertexData& dt, const Point& p);

//Finds edge by two vertices.
class EdgeFinder{
	vector<Connectivity::VertexEdgeR> ve;
	const EdgeData* data;
public:
	EdgeFinder(const EdgeData& data);

	//returns:
	// <0> found edge or nullptr
	// <1> correct direction of [v1, v2]
	std::tuple<shared_ptr<Edge>, bool>              
	find(Vertex* v1, Vertex* v2);
};

//vertex match searcher.
//ids features are not touched by this implementation
class VertexMatch{
	VertexData srt;
public:
	VertexMatch(const VertexData& vd);
	shared_ptr<Vertex> find(const Point& p);
	VertexData find(const vector<Point>& p);
};

class RasterizeEdges{
	std::shared_ptr<BoundingBoxFinder> bbf;
	//true for squares which contain an edge
	vector<bool> black_squares;
public:
	RasterizeEdges(const EdgeData& ed, const BoundingBox& bb, double step);

	// 0 for black squares
	// 1, 2, 3 for grouped white squares
	// if use_groups = false, all whites have unity feature
	vector<int> colour_squares(bool use_groups) const;

	BoundingBoxFinder& bbfinder(){ return *bbf; }
};

}

namespace Contour{ namespace Finder{

//->(INSIDE, BOUND, OUTSIDE). Direction is ignored.
//bb is a helper bounding box of input contour. 
int WhereIs(const EdgeData&, const Point& p, BoundingBox* bb=0);

// ======== Crosses and intersections
//finds first cross (with respect to length of c1) of contours c1, c2.
//returns <0>: if cross was found
//        <1>: cross point
//        <2,3>: normalized length coordinate of intersection
std::tuple<bool, Point, double, double> 
Cross(const EdgeData& c1, const EdgeData& c2);

vector<std::tuple<bool, Point, double, double>>
CrossAll(const EdgeData& c1, const EdgeData& c2);

//finds first self cross (with respect to length of c1)
//returns <0>: if cross was found
//        <1>: cross point
//        <2,3>: normalized section length coordinates of intersection
//        <4,5>: local indices of crossed sections
std::tuple<bool, Point, double, double, int, int>
SelfCross(const EdgeData& c1);

//Calculate points position with respect to contour.
//Contour direction is taken into account.
//Builds vector with INSIDE/OUTSIDE/BOUND for each point.
vector<int> SortOutPoints(const EdgeData& t1, const vector<Point>& pnt);
vector<int> SortOutPoints(const Contour::Tree& t1, const vector<Point>& pnt);

}}

}

#endif
