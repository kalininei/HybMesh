#ifndef  SURFACE_GRID3D_HPP
#define  SURFACE_GRID3D_HPP

#include "primitives_grid3d.hpp"

namespace HMGrid3D{
struct SurfaceTree;

struct Surface{
	//data
	ShpVector<Face> faces;

	//features
	ShpVector<Face> allfaces() const;
	ShpVector<Edge> alledges() const;
	ShpVector<Vertex> allvertices() const;
	int n_faces() const { return faces.size(); }

	// ================= Builders
	//reversetp = 
	//  1: reverse faces so that all grid cells be on their right sides.
	//  0: do not reverse
	// -1: reverse faces so that all grid cells be on their left sides.
	// makes a shallow copy of grid data, so reversing will alter grid structure
	static Surface FromBoundaryType(HMGrid3D::GridData& g, int btype, int reversetp);
	static std::map<int, shared_ptr<Surface>> FromBoundaryType(HMGrid3D::GridData& g, int reversetp);

	static SurfaceTree BuildTree(HMGrid3D::GridData& g, int reversetp);

	//extracts subsurface which contains v
	static Surface SubSurface(const Surface& s, Vertex* v);

	//topologically unique rearrange with respect to given edge
	//edges and faces data would be reordered, face vector will be permuted.
	static void FaceRearrange(Surface& s, const Edge* ed);

	// ================= Algos
	static bool MatchTopology(const Surface& a, const Surface& b);

	//!!! all boundary edges must have correct direction
	static ShpVector<Edge> ExtractBoundary(const Surface& a, Vertex v);
};

struct SurfaceTree{
	struct TNode: public Surface{
		TNode():level(0){}
		TNode(const Surface& s): Surface(s), level(0){}
		weak_ptr<TNode> parent;
		WpVector<TNode> children;
		int level;
	};
	ShpVector<TNode> nodes;

	//reallocates nodes. all data remains.
	SurfaceTree reallocate_nodes() const;

	vector<SurfaceTree> crop_level1() const;
	static SurfaceTree Substract(const SurfaceTree& from, const SurfaceTree& what);
	static SurfaceTree Assemble(const vector<Surface>& data);
};

}

#endif

