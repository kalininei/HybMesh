#ifndef HMCONT2D_TREE_HPP
#define HMCONT2D_TREE_HPP

#include "primitives2d.hpp"
#include "contour.hpp"

namespace HM2D{ namespace Contour{

struct Tree{
	Tree(const vector<EdgeData>& dt={}) { for (auto& it: dt) add_contour(it); }
	struct TNode{
		EdgeData contour;
		weak_ptr<TNode> parent;
		WpVector<TNode> children;
		int level;

		TNode():level(0){}
		TNode(const EdgeData& s, int level=0): contour(s), level(level){}
		TNode(EdgeData&& s, int level=0): contour(std::move(s)), level(level){}

		bool isdetached()const{ return level<0; }
		bool isbound()const{ return level>=0; }
		bool isroot()const{ return level==0; }
		bool isouter()const{ return level % 2 == 0;}
		bool isinner()const{ return level % 2 == 1;}
	};
	ShpVector<TNode> nodes;

	//features
	ShpVector<TNode> roots() const;
	ShpVector<TNode> detached_contours() const;
	ShpVector<TNode> bound_contours() const;
	EdgeData alledges() const;
	EdgeData alledges_bound() const;
	shared_ptr<TNode> find_node(const Point* p) const;
	shared_ptr<TNode> find_node(const Edge* p) const;

	//builders
	void add_contour(const EdgeData& cont);
	void add_detached_contour(const EdgeData& cont);
	void add_contour(EdgeData&& cont);
	void add_detached_contour(EdgeData&& cont);
	void remove_contour(TNode* n);

	//using node->parents calculate children and nesting level
	void update_topology();
	void remove_detached();

	//functions
	double area() const;
	int whereis(const Point&) const;

	//assemble procedure: returned tree uses edges of data
	static Tree Assemble(const EdgeData& data);
	static Tree GridBoundary(const GridData& data);
	//first node of each returning entry is a single root contour.
	//all others are its children
	static vector<Tree> GridBoundary01(const GridData& data);

	//returns vector of trees each of which contains single root node
	//and its children.
	static vector<Tree> CropLevel01(const Tree& data);

	//lv = 3 => deepcopy all primitives
	//lv = 2 => leave vertices (realloc tnodes and edges)
	//lv = 1 => leave vertices and edges (realloc tnodes)
	//lv = 0 => leave everything: vertices, edges and tnodes
	static Tree DeepCopy(const Tree& from, int level=3);

};


}}




#endif
