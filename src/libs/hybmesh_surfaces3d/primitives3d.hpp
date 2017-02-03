#ifndef PRIMITIVES_GRID3D_HPP
#define PRIMITIVES_GRID3D_HPP
#include "hmproject.h"
#include "bgeom2d.h"
#include "bgeom3d.h"

namespace HM3D{

struct Vertex;
struct Edge;
struct Face;
struct Cell;
typedef ShpVector<Vertex> VertexData;
typedef ShpVector<Edge> EdgeData;
typedef ShpVector<Face> FaceData;
typedef ShpVector<Cell> CellData;


struct Vertex: public Point3{
	//==== data
	mutable int id;

	//==== constructor
	Vertex(double x=0, double y=0, double z=0): Point3(x, y, z){}
	Vertex(const Point3& p): Point3(p){}
};

struct Edge{
	// ==== Data
	std::array<shared_ptr<Vertex>, 2> vertices;
	mutable int id;
	
	//==== Constructor
	Edge(shared_ptr<Vertex> p1, shared_ptr<Vertex> p2): vertices {p1, p2}{}
	Edge(){}

	//==== Features
	shared_ptr<Vertex> first() const { return vertices[0]; }
	shared_ptr<Vertex> last() const { return vertices[1]; }
	double measure() const;
	double length() const;

	//==== Funcs
	void reverse();

};

struct Face{
	// ===== Data
	//edges are sorted
	EdgeData edges;
	int boundary_type;
	//if one looks at the face and sees it in a couterclockwise order 
	//then he looks from a right direction.
	std::weak_ptr<Cell> left, right;
	mutable int id;

	// ===== Constructor
	Face(const ShpVector<Edge>& e={}, int bt=0): edges(e), boundary_type(bt){};

	// ===== Features
	bool is_boundary() const { return left.expired() || right.expired(); }
	bool has_right_cell() const { return !right.expired(); }
	bool has_left_cell() const { return !left.expired(); }
	std::array<Point3, 3> mean_points() const;
	Vect3 left_normal() const;
	//first vertex is common to last and first edge
	ShpVector<Vertex> sorted_vertices() const;

	// ===== funcs
	void reverse();
	void correct_edge_directions();
	bool is_positive_edge(int eindex);

	// ==== change data
	// changes edge entry in this->edges with underlying vertices
	// if from was not found in edges returns false
	// from and to edges should be constructed by equal valued end vertices
	bool change_edge(shared_ptr<Edge> from, shared_ptr<Edge> to);

};

struct Cell{
	// ===== Data
	FaceData faces;
	mutable int id;

	// ==== Features
	//number of faces, edges, vertices
	std::tuple<int, int, int> n_fev() const;
	//cell volume
	double volume() const;

	// ==== change data
	// changes face entry in this->faces with underlying edges and vertices
	// doesn't change input faces left/right attibutes,
	// if from was not found in faces returns false
	// from and to faces should be constructed by equal valued vertices
	// and have same dirction
	bool change_face(shared_ptr<Face> from, shared_ptr<Face> to);
};

//Grid primitives collection
struct GridData{
	VertexData vvert;
	EdgeData vedges;
	FaceData vfaces;
	CellData vcells;

	// ====== methods
	void clear(){
		vvert.clear();
		vedges.clear();
		vfaces.clear();
		vcells.clear();
	}
	//set id's of primitives to its actual indices
	void enumerate_all() const;
};

//deep copy procedures
void DeepCopy(const VertexData& from, VertexData& to);
void DeepCopy(const EdgeData& from, EdgeData& to, int level=1);
void DeepCopy(const FaceData& from, FaceData& to, int level=2);
void DeepCopy(const CellData& from, CellData& to, int level=3);
void DeepCopy(const GridData& from, GridData& to, int level=3);

//extract procedures
VertexData AllVertices(const EdgeData& from);
VertexData AllVertices(const FaceData& from);
VertexData AllVertices(const CellData& from);
EdgeData AllEdges(const FaceData& from);
EdgeData AllEdges(const CellData& from);
FaceData AllFaces(const CellData& from);
std::tuple<VertexData> AllPrimitives(const EdgeData& from);
std::tuple<VertexData, EdgeData> AllPrimitives(const FaceData& from);
std::tuple<VertexData, EdgeData, FaceData> AllPrimitives(const CellData& from);

//miscellaneous algorithms
std::tuple<
	Vertex*,     //pointer to closest vertex
	int,         //index of closest vertex within vec
	double       //measure to closest vertex
> FindClosestVertex(const VertexData& vec, Vertex v);

//connect edges at data starting from vertex closest to v till close or end reached.
//chooses direction according to edge which includes v as start vertex
//!!! If edge vertex has more then 2 connections or 
//    multiple edges include v as start point result is undefined (false assert in debug mode)
//!!! Doesn't go backward
EdgeData Connect(const EdgeData& data, Vertex v);

//cells volumes
vector<double> Volumes(const CellData& cd);
double SumVolumes(const CellData& cd);

//split procedures
vector<FaceData> SplitData(const FaceData& data);
vector<CellData> SplitData(const CellData& data);
vector<GridData> SplitData(const GridData& data);

//scaling procedures
ScaleBase3 Scale01(VertexData&, double a=1);
ScaleBase3 Scale01(FaceData&, double a=1);
void Scale(VertexData&, const ScaleBase3& sc);
void Scale(FaceData&, const ScaleBase3& sc);
void Unscale(VertexData&, const ScaleBase3& sc);
void Unscale(FaceData&, const ScaleBase3& sc);

}

#endif

