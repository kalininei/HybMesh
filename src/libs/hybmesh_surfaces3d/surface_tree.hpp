#ifndef HYBMESH_SURFACE_TREE_HPP
#define HYBMESH_SURFACE_TREE_HPP
#include "surface.hpp"

namespace HM3D{ namespace Surface{

struct Tree{
	struct TNode{
		FaceData surface;
		weak_ptr<TNode> parent;
		WpVector<TNode> children;
		int level;

		TNode():level(0){}
		TNode(const FaceData& s, int level=0): surface(s), level(level){}
		bool isopen(){ return level<0; }
		bool isclosed(){ return level>=0; }
	};
	ShpVector<TNode> nodes;
	ShpVector<TNode> roots() const;
	FaceData allfaces() const;

	vector<Tree> crop_level1() const;
	static Tree Assemble(const FaceData& data);
};


}}

#endif
