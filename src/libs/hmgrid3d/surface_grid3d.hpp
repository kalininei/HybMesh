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

	// ================= Builders
	//reversetp = 
	//  1: reverse faces so that all grid cells be on their right sides.
	//  0: do not reverse
	// -1: reverse faces so that all grid cells be on their left sides.
	// makes a shallow copy of grid data, so reversing will alter grid structure
	static Surface FromBoundaryType(HMGrid3D::GridData& g, int btype, int reversetp);
	static std::map<int, Surface> ByBoundaryTypes(HMGrid3D::GridData& g, int reversetp);
	static std::map<int, Surface> ByBoundaryTypes(const HMGrid3D::GridData& g);
	static Surface GridSurface(HMGrid3D::GridData& g, int reversetp);
	static Surface GridSurface(const HMGrid3D::GridData& g);

	//extracts subsurface which contains v
	static Surface SubSurface(const Surface& s, Vertex* v);
	//all subsurfaces
	static vector<Surface> AllSubSurfaces(const Surface& s);
	//extract smooth surface sections with normal deviations less than angle(deg)
	static vector<Surface> ExtractSmooth(const Surface& s, double angle);

	//topologically unique rearrange with respect to given edge
	//edges and faces data would be reordered, face vector will be permuted.
	static void FaceRearrange(Surface& s, const Edge* ed);

	// ================= Algos
	static bool MatchTopology(const Surface& a, const Surface& b);

	//!!! all boundary edges should have correct direction
	static ShpVector<Edge> ExtractBoundary(const Surface& a, Vertex v);
	//return[0] - positive closed boundaries, return[1] - negative
	//!!! face directions should match in surface
	static std::array<vector<EdgeData>, 2>
		ExtractAllBoundaries(const Surface& a, Vect3 right_normal);

	//signed volume starting from first three points of first face
	//!!! surface should bound singly connected domain
	//!!! all faces should have same direction relative to volume
	static double Volume(const FaceData& fc);
	static double Volume(const Surface& a){ return Volume(a.faces); }
};

struct SurfaceTree{
	struct TNode: public Surface{
		TNode():level(0){}
		TNode(const Surface& s, int level=0): Surface(s), level(level){}
		weak_ptr<TNode> parent;
		WpVector<TNode> children;
		int level;
		bool isopen(){ return level<0; }
		bool isclosed(){ return level>=0; }
	};
	ShpVector<TNode> nodes;
	ShpVector<TNode> roots() const;

	HMGrid3D::EdgeData alledges() const;
	HMGrid3D::VertexData allvertices() const;
	HMGrid3D::FaceData allfaces() const;

	vector<SurfaceTree> crop_level1() const;
	static SurfaceTree Assemble(const Surface& data);
};

//temporary reverts all face edges direction so that all
//face dirctions match each other.
//When this object dies it reverts faces back to its original state.
struct SurfTReverterBase{
	bool reverted;
	SurfTReverterBase(): reverted(false){}
	virtual ~SurfTReverterBase(){}
	void revert(){
		if (!reverted) _do_revert();
		reverted=true;
	}
	void revert_back(){
		if (reverted) _undo_revert();
		reverted=false;
	}
private:
	virtual void _do_revert() = 0;
	virtual void _undo_revert() = 0;
};
struct SurfTReverter: public SurfTReverterBase{
	Surface* obj;
	vector<bool> need_revert;
	//using const since all surf changes are temporal
	SurfTReverter(const Surface& srf);
	~SurfTReverter(){ revert_back(); }

	void reverse_all();
private:
	void _do_revert() override;
	void _undo_revert() override;
};

//tree internal area is located to the left of even leveled surfaces
//and to the right of odd leveled surfaces.
struct SurfTreeTReverter: public SurfTReverterBase{
	SurfaceTree* obj;
	ShpVector<SurfTReverter> openrevs;
	ShpVector<SurfTReverter> closedrevs;
	//using const since all surface changes are temporal
	SurfTreeTReverter(const SurfaceTree& srf);
	~SurfTreeTReverter(){ revert_back(); }
private:
	void _do_revert() override;
	void _undo_revert() override;
};

}

#endif

