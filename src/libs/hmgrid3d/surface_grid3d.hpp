#ifndef SURFACE_GRID3D_HPP
#define SURFACE_GRID3D_HPP

#include "primitives_grid3d.hpp"

namespace HMGrid3D{

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
	static std::map<int, Surface> ByBoundaryTypes(HMGrid3D::GridData& g, int reversetp);
	static std::map<int, Surface> ByBoundaryTypes(const HMGrid3D::GridData& g);

	//extracts subsurface which contains v
	static Surface SubSurface(const Surface& s, Vertex* v);

	//topologically unique rearrange with respect to given edge
	//edges and faces data would be reordered, face vector will be permuted.
	static void FaceRearrange(Surface& s, const Edge* ed);

	// ================= Algos
	static bool MatchTopology(const Surface& a, const Surface& b);

	//!!! all boundary edges should have correct direction
	static ShpVector<Edge> ExtractBoundary(const Surface& a, Vertex v);
};

}

#endif

