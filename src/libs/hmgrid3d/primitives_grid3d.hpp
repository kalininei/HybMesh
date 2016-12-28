#ifndef PRIMITIVES_GRID3D_HPP
#define PRIMITIVES_GRID3D_HPP
#include "hmproject.h"
#include "bgeom2d.h"
#include "bgeom3d.h"

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
//class which keeps enumeration of *Data
//and places it back on delete
template<class C>
struct RestoreIds{
	vector<int> ids;
	vector<typename C::value_type> _deepdata;
	const C* _shallowdata;

	//use deepcopy=true only if entries of pvec
	//could be reallocated/destructed between constructing and destructing of this object
	//In the latter case this object restores ids of old pvec objects anyway
	RestoreIds(const C& pvec, bool deepcopy=false): _shallowdata(&pvec){
		ids.resize(pvec.size());
		for (int i=0; i<pvec.size(); ++i) ids[i] = pvec[i]->id;
		if (deepcopy){
			_deepdata.resize(pvec.size());
			for (int i=0; i<pvec.size(); ++i) _deepdata[i] = pvec[i];
		}
	}
	~RestoreIds(){
		if (_deepdata.size() == 0)
			for (int i=0; i<ids.size(); ++i) (*_shallowdata)[i]->id = ids[i];
		else
			for (int i=0; i<ids.size(); ++i) _deepdata[i]->id = ids[i];
	}
};

struct Vertex: public Point3{
	//==== data
	mutable int id;

	//==== constructor
	Vertex(double x=0, double y=0, double z=0): Point3(x, y, z){}
	Vertex(const Point3& p): Point3(p){}

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
		struct EndVertexEdgeR{
			shared_ptr<Vertex> v;
			vector<int> eind;
			size_t size() const { return eind.size(); }
		};
		static vector<EndVertexEdgeR> EndVertexEdge(const ShpVector<Edge>& data);
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
	//first vertex is common to last and first edge
	ShpVector<Vertex> sorted_vertices() const;
	ShpVector<Edge> alledges() const;
	ShpVector<Vertex> allvertices() const;

	// ===== Constructor
	Face(const ShpVector<Edge>& e={}, int bt=0): edges(e), boundary_type(bt){};

	// ===== funcs
	void reverse();
	void correct_edge_directions();
	bool is_positive_edge(int eindex);
	std::array<Point3, 3> mean_points() const;
	Vect3 left_normal() const;

	// ==== change data
	// changes edge entry in this->edges with underlying vertices
	// if from was not found in edges returns false
	// from and to edges should be constructed by equal valued end vertices
	bool change_edge(shared_ptr<Edge> from, shared_ptr<Edge> to);

	// ===== Algos
	static std::vector<ShpVector<Face>> SubDivide(const ShpVector<Face>& fvec);
	//assign boundary by centeral vertex and old boundary type.
	static void SetBoundaryTypes(const ShpVector<Face>& fvec, std::function<int(Vertex, int)> bfun);

	// ===== Connectivity tables
	struct Connectivity{
		struct EdgeFaceR{
			shared_ptr<Edge> e;
			vector<int> find;  //faces indicies
			size_t size() const { return find.size(); }
		};
		static vector<EdgeFaceR> EdgeFace(const FaceData& data);

		//for each edge- > find - face index,
		//                 locind - local edge index within the face,
		//                 posdir - is edge directed according to face
		struct EdgeFaceExtendedR: public EdgeFaceR{
			vector<int> locind;
			vector<bool> posdir;
		};
		static vector<EdgeFaceExtendedR> EdgeFaceExtended(const FaceData& data);

		static vector<vector<int>>
		FaceFace(const FaceData& data);
		static vector<vector<int>>
		FaceFace(const vector<EdgeFaceR>& edge_face, int nfaces);
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

	// ==== change data
	// changes face entry in this->faces with underlying edges and vertices
	// doesn't change input faces left/right attibutes,
	// if from was not found in faces returns false
	// from and to faces should be constructed by equal valued vertices
	// and have same dirction
	bool change_face(shared_ptr<Face> from, shared_ptr<Face> to);
	
	// ==== additional information
	double volume() const;
	static vector<double> Volumes(const CellData& cd);
	static double SumVolumes(const CellData& cd);
	
	// ===== Data access
	ShpVector<Vertex> allvertices() const;
	ShpVector<Face> allfaces() const;
	ShpVector<Edge> alledges() const;
};

//deep copy procedures
void DeepCopy(const VertexData& from, VertexData& to);
void DeepCopy(const EdgeData& from, EdgeData& to, int level=1);
void DeepCopy(const FaceData& from, FaceData& to, int level=2);
void DeepCopy(const CellData& from, CellData& to, int level=3);
//extract procedures
VertexData AllVertices(const EdgeData& from);
VertexData AllVertices(const FaceData& from);
VertexData AllVertices(const CellData& from);
VertexData AllEndVertices(const EdgeData& from);
VertexData AllEndVertices(const FaceData& from);
VertexData AllEndVertices(const CellData& from);

EdgeData AllEdges(const FaceData& from);
EdgeData AllEdges(const CellData& from);

FaceData AllFaces(const CellData& from);

std::tuple<VertexData> AllPrimitives(const EdgeData& from);
std::tuple<VertexData, EdgeData> AllPrimitives(const FaceData& from);
std::tuple<VertexData, EdgeData, FaceData> AllPrimitives(const CellData& from);


//Grid primitives collection
struct GridData{
	VertexData vvert;
	EdgeData vedges;
	FaceData vfaces;
	CellData vcells;

	// ====== methods
	void clear(){ vvert.clear(); vedges.clear(); vfaces.clear(); vcells.clear(); }
	//set id's of primitives to its actual indicies
	void enumerate_all() const;

};



}

#endif

