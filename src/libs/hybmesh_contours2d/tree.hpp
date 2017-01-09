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
		bool isdetached(){ return level<0; }
		bool isbound(){ return level>=0; }
	};
	ShpVector<TNode> nodes;

	//features
	ShpVector<TNode> roots() const;
	ShpVector<TNode> detached_contours() const;
	ShpVector<TNode> bound_contours() const;
	EdgeData alledges() const;
	shared_ptr<TNode> find_node(const Point* p) const;
	shared_ptr<TNode> find_node(const Edge* p) const;

	//builders
	void add_contour(const EdgeData& cont);
	void add_detached_contour(const EdgeData& cont);
	//using node->parents calculate children and nesting level
	void update_topology();
	void remove_detached();

	//functions
	double area() const;
	int whereis(const Point&) const;

	//assemble procedure: returned tree uses edges of data
	static Tree Assemble(const EdgeData& data);

	//lv = 3 => deepcopy all primitives
	//lv = 2 => leave vertices (realloc tnodes and edges)
	//lv = 1 => leave vertices and edges (realloc tnodes)
	//lv = 0 => leave everything: vertices, edges and tnodes
	static Tree DeepCopy(const Tree& from, int level=3);
};


}}




#endif
