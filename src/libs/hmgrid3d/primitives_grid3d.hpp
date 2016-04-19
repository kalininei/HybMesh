#ifndef PRIMITIVES_GRID3D_HPP
#define PRIMITIVES_GRID3D_HPP
#include "hmproject.h"
#include "bgeom2d.h"

namespace HMGrid3D{

struct Vertex;
struct Edge;
struct Face;
struct Cell;
typedef ShpVector<Vertex> VertexData;
typedef ShpVector<Edge> EdgeData;
typedef ShpVector<Face> FaceData;
typedef ShpVector<Cell> CellData;

template<class C>
void enumerate_ids_pvec(const C& inp){
	for (int i=0; i<inp.size(); ++i) inp[i]->id = i;
}
template<class C>
void constant_ids_pvec(const C& inp, int val){
	for (int i=0; i<inp.size(); ++i) inp[i]->id = val;
}

struct Vertex{
	//==== data
	double x, y, z;
	mutable int id;

	//==== constructor
	Vertex(double _x=0, double _y=0, double _z=0): x(_x), y(_y), z(_z){}

	//==== features
	static double measure(const Vertex& a, const Vertex& b);

	//===== algos
	static std::tuple<
		Vertex*,     //pointer to closest vertex
		int,         //index of closest vertex within vec
		double       //measure to closest vertex
	> FindClosestVertex(const ShpVector<Vertex>& vec, Vertex v);

	static void Unscale2D(ShpVector<Vertex>& vec, ScaleBase sc);
};

struct Edge{
	// ==== Data
	ShpVector<Vertex> vertices;
	mutable int id;
	
	//==== Constructor
	Edge(shared_ptr<Vertex> p1, shared_ptr<Vertex> p2): vertices {p1, p2}{}
	Edge(){}

	//==== Features
	shared_ptr<Vertex> first() const { return vertices[0]; }
	shared_ptr<Vertex> last() const { return vertices.back(); }
	double measure() const;
	double length() const;

	//==== Funcs
	void reverse();

	//==== Algos
	//connect edges at data starting from vertex closest to v till close or end reached.
	//chooses direction according to edge which includes v as start vertex
	//!!! If edge vertex has more then 2 connections or 
	//    multiple edges include v as start point result is undefined (false assert in debug mode)
	//!!! Doesn't go backward
	static ShpVector<Edge> Connect(const ShpVector<Edge>& data, Vertex v);

	// ===== Connectivity tables
	struct Connectivity{
		static std::map<shared_ptr<Vertex>, vector<int>>
		EndVertexEdge(const ShpVector<Edge>& data);
	};
};

struct CurvEdge: public Edge{
	double curv;
	CurvEdge(shared_ptr<Vertex> p1, shared_ptr<Vertex> p2, double c): Edge(p1, p2), curv(c){}
	CurvEdge(){}
};

struct Face{
	// ===== Data
	//edges are sorted
	ShpVector<Edge> edges;
	int boundary_type;
	//if one looks at the face and sees it in a couterclockwise order 
	//then he looks from a right direction.
	std::weak_ptr<Cell> left, right;
	mutable int id;

	// ===== Features
	bool is_boundary() const { return left.expired() || right.expired(); }
	bool has_right_cell() const { return !right.expired(); }
	bool has_left_cell() const { return !left.expired(); }
	int n_edges() const { return edges.size(); }

	// ===== Data access
	ShpVector<Vertex> sorted_vertices() const;
	ShpVector<Edge> alledges() const;
	ShpVector<Vertex> allvertices() const;

	// ===== Constructor
	Face(const ShpVector<Edge>& e={}, int bt=0): edges(e), boundary_type(bt){};

	// ===== funcs
	void reverse();
	void correct_edge_directions();

	// ===== Algos
	static std::vector<ShpVector<Face>> SubDivide(const ShpVector<Face>& fvec);
	//assign boundary by centeral vertex and old boundary type.
	static void SetBoundaryTypes(const ShpVector<Face>& fvec, std::function<int(Vertex, int)> bfun);

	// ===== Connectivity tables
	struct Connectivity{
		static std::map<shared_ptr<Edge>, vector<int>>
		EdgeFace(const ShpVector<Face>& data);
		static vector<vector<int>>
		FaceFace(const ShpVector<Face>& data);
		static vector<vector<int>>
		FaceFace(const std::map<shared_ptr<Edge>, vector<int>>& edge_face, int nfaces);
	};
};

struct Cell{
	// ===== Data
	ShpVector<Face> faces;
	mutable int id;

	// ==== Features
	int n_faces() const;
	int n_edges() const;
	int n_vertices() const;
	std::tuple<int, int, int> n_fev() const; //number of faces, edges, vertices

	// ===== Data access
	ShpVector<Vertex> allvertices() const;
	ShpVector<Face> allfaces() const;
	ShpVector<Edge> alledges() const;
};

//Grid primitives collection
struct GridData{
	VertexData vvert;
	EdgeData vedges;
	FaceData vfaces;
	CellData vcells;

	// ====== methods
	//set id's of primitives to its actual indicies
	void enumerate_all() const;
};



}

#endif

