#ifndef HYBMESH_PRIMITIVES2D_HPP
#define HYBMESH_PRIMITIVES2D_HPP

#include "hmproject.h"
#include "bgeom2d.h"

namespace HM2D{

struct Vertex;
struct Edge;
struct Cell;
typedef ShpVector<Vertex> VertexData;
typedef ShpVector<Edge> EdgeData;
typedef ShpVector<Cell> CellData;

struct Vertex: public Point{
	//==== data
	mutable int id;

	//==== constructor
	Vertex(double x=0, double y=0): Point(x, y){}
	Vertex(const Point& p): Point(p){}
};

struct Edge{
	// ==== Data
	std::array<shared_ptr<Vertex>, 2> vertices;
	mutable int id;
	int boundary_type;
	weak_ptr<Cell> left, right;
	
	//==== Constructor
	Edge(shared_ptr<Vertex> p1, shared_ptr<Vertex> p2): vertices {p1, p2}, boundary_type(0){}
	Edge(): boundary_type(0){}

	//==== Features
	shared_ptr<Vertex> first() const { return vertices[0]; }
	shared_ptr<Vertex> last() const { return vertices[1]; }
	bool has_right_cell() const { return !right.expired(); }
	bool has_left_cell() const { return !left.expired(); }
	bool is_boundary() const { return left.expired() || right.expired(); }
	double measure() const { return Point::meas(*vertices[0], *vertices[1]); }
	double length() const { return sqrt(measure()); }
	bool connected_to(const Edge& e2) const{
		return first() == e2.last() ||
		       first() == e2.first() ||
		       last() == e2.first() ||
		       last() == e2.last();
	}
	shared_ptr<Vertex> sibling(const Point* p1) const{
		if (p1 == vertices[0].get()) return vertices[1];
		else return vertices[0];
	}
	Point center() const { return Point::Weigh(*vertices[0], *vertices[1], 0.5); }

	//==== Funcs
	void reverse(){
		std::swap(vertices[0], vertices[1]);
		std::swap(left, right);
	}
};

struct Cell{
	// ===== Data
	//edges are sorted in counterclockwise direction
	EdgeData edges;
	mutable int id;

	// ===== Constructor
	Cell(const ShpVector<Edge>& e={}): edges(e){};

	// ===== Features
	//first vertex is common to last and first edge
	ShpVector<Vertex> sorted_vertices() const;
	//signed area. should be positive for valid cells
	double area() const;
	//number of edges, vertices
	std::tuple<int, int> n_ev() const;

	// ===== funcs
	void correct_edge_directions();
	bool is_positive_edge(int eindex);
};

//Grid primitives collection
struct GridData{
	VertexData vvert;
	EdgeData vedges;
	CellData vcells;

	// ====== methods
	void clear(){
		vvert.clear();
		vedges.clear();
		vcells.clear();
	}
	//set id's of primitives to its actual indices
	void enumerate_all() const;
};

//deep copy procedures
void DeepCopy(const VertexData& from, VertexData& to);
void DeepCopy(const EdgeData& from, EdgeData& to, int level=1);
void DeepCopy(const CellData& from, CellData& to, int level=2);
void DeepCopy(const GridData& from, GridData& to, int level=2);

//contains procedures: returns data object or nullptr
shared_ptr<Vertex> Contains(const VertexData& data, const Point* pnt);
shared_ptr<Vertex> Contains(const EdgeData& data, const Point* pnt);
shared_ptr<Edge> Contains(const EdgeData& data, const Edge* ed);
shared_ptr<Vertex> Contains(const CellData& data, const Point* pnt);
shared_ptr<Edge> Contains(const CellData& data, const Edge* ed);
shared_ptr<Cell> Contains(const CellData& data, const Cell* c);

//extract procedures
VertexData AllVertices(const EdgeData& from);
VertexData AllVertices(const CellData& from);
EdgeData AllEdges(const CellData& from);
std::tuple<VertexData> AllPrimitives(const EdgeData& from);
std::tuple<VertexData, EdgeData> AllPrimitives(const CellData& from);

double Length(const EdgeData&);
vector<double> ELengths(const EdgeData&);

//closest edge-> returns 
//<0> edge index, -1 if dt is empty
//<1> distance,
//<2> weight of closest point within edge.
std::tuple<int, double, double>
FindClosestEdge(const EdgeData& dt, const Point& p);

//closest point lying on edge
Point FindClosestEPoint(const EdgeData& dt, const Point& p);

//pointer to closest point:
//<0> - point index
//<1> - squared distance to point
std::tuple<int, double> FindClosestNode(const VertexData& dt, const Point& p);


BoundingBox BBox(const VertexData&, double eps=0.);
BoundingBox BBox(const EdgeData&, double eps=0.);
BoundingBox BBox(const CellData&, double eps=0.);

ScaleBase Scale01(EdgeData&);
ScaleBase Scale01(VertexData&);
void Scale(EdgeData&, const ScaleBase& sc);
void Unscale(EdgeData&, const ScaleBase& sc);

}

#endif
