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

		bool isdetached()const{ return level<0; }
		bool isbound()const{ return level>=0; }
		bool isroot()const{ return level==0; }
		bool isouter()const{ return level % 2 == 0;}
		bool isinner()const{ return level % 2 == 1;}
	};
	ShpVector<TNode> nodes;
	ShpVector<TNode> roots() const;
	FaceData allfaces() const;

	static vector<Tree> CropLevel01(const Tree&);
	static Tree Assemble(const FaceData& data);
};


}}

#endif
