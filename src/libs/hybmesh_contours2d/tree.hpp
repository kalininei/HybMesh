#ifndef HMCONT2D_TREE_HPP
#define HMCONT2D_TREE_HPP

#include "primitives2d.hpp"
#include "contour.hpp"

namespace HM2D{ namespace Contour{

struct Tree{
	struct TNode{
		EdgeData contour;
		weak_ptr<TNode> parent;
		WpVector<TNode> children;
		int level;

		TNode():level(0){}
		TNode(const EdgeData& s, int level=0): contour(s), level(level){}
		bool isopen(){ return level<0; }
		bool isclosed(){ return level>=0; }
	};
	ShpVector<TNode> nodes;
	ShpVector<TNode> roots() const;
	ShpVector<TNode> open_contours() const;
	ShpVector<TNode> closed_contours() const;
	//features
	EdgeData alledges() const;
	double area() const;
	int whereis(const Point&) const;
	shared_ptr<TNode> find_node(const Point* p) const;
	shared_ptr<TNode> find_node(const Edge* p) const;

	//builders
	void AddContour(const EdgeData& cont);
	//using node->parents calculate children + checks direction
	void UpdateTopology();

	//procedures
	vector<Tree> crop_level1() const;
	void remove_opens();
	static Tree Assemble(const EdgeData& data);

	//lv = 3 => deepcopy all primitives
	//lv = 2 => leave vertices (realloc tnodes and edges)
	//lv = 1 => leave vertices and edges (realloc tnodes)
	//lv = 0 => leave vertices, edges and tnodes
	static Tree DeepCopy(const Tree& from, int level=3);
};


}}




#endif
