#ifndef PRIMITIVES_GRID3D_HPP
#define PRIMITIVES_GRID3D_HPP
#include "hmproject.h"

namespace HMGrid3D{

struct Vertex;
struct Edge;
struct Face;
struct Cell;
struct Grid;

struct Vertex{
	double x, y, z;
	Vertex(double _x=0, double _y=0, double _z=0): x(_x), y(_y), z(_z){}
};

struct Edge{
	// ==== Data
	ShpVector<Vertex> vertices;
	
	// ==== constructor
	Edge(shared_ptr<Vertex> p1, shared_ptr<Vertex> p2): vertices {p1, p2}{}

	// ==== features
	double measure() const;
	double length() const;
};

struct Face{
	// ===== Data
	//edges are sorted
	ShpVector<Edge> edges;
	int boundary_type;
	//if one looks at the face and sees it in a couterclockwise order 
	//then he looks from a right direction.
	shared_ptr<Cell> left, right;

	// ===== Features
	bool is_boundary() const { return left==0 || right==0; }

	// ===== Data access
	ShpVector<Vertex> sorted_vertices() const;
	ShpVector<Vertex> allvertices() const;

	// ===== Constructor
	Face(const ShpVector<Edge>& e={}, int bt=0): edges(e), boundary_type(bt){};
};

struct Cell{
	// ===== Data
	ShpVector<Face> faces;

	// ===== Data access
	ShpVector<Vertex> allvertices() const;
	ShpVector<Face> allfaces() const;
	ShpVector<Edge> alledges() const;
};

struct Grid{
	// ====== Data
	ShpVector<Cell> cells;

	// ===== Features
	int n_cells() const;
	int n_faces() const;
	int n_edges() const;
	int n_vertices() const;
	// ===== Data access
	ShpVector<Vertex> allvertices() const;
	ShpVector<Edge> alledges() const;
	ShpVector<Face> allfaces() const;
	ShpVector<Cell> allcells() const;
};

}

#endif

