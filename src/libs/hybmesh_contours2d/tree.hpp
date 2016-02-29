#ifndef HMCONT2D_TREE_HPP
#define HMCONT2D_TREE_HPP

#include "collections.hpp"
#include "contour.hpp"

namespace HMCont2D{

// ============== Tree: set of non-overlapping contours
struct ContourTree: public ECollection {
	//All contours should be closed
	//edges are shared between inherited this->data and
	//contours[i].data;

	//contour with its relations to others
	struct TreeNode: public Contour{
		TreeNode* parent;
		vector<TreeNode*> children;
	};

	//all nodes
	ShpVector<TreeNode> nodes;

	//get roots node: always outer contours
	vector<TreeNode*> roots() const;
	//number of contours
	virtual int cont_count() const { return nodes.size(); }

	//get contour by point
	virtual Contour* get_contour(Point* p) const;
	//get contour by edge
	virtual Contour* get_contour(Edge* p) const;
	//get contour by index
	virtual Contour* get_contour(int i) const;

	//Methods
	virtual void AddContour(shared_ptr<Contour>& c);
	void AddContour(const Contour& c);
	virtual void RemovePoints(const vector<const Point*>& p);

	//Returns true if point lies strictly within/without closed contour
	//for points on edges result is unpredirctable
	bool IsWithin(const Point& p) const;
	bool IsWithout(const Point& p) const;

	//using node->parents calculate children + checks direction
	void UpdateTopology();
	//fills this->data with data from nodes
	//should be called if nodes edges are added or deleted to make this->data actual
	virtual void ReloadEdges();

	//Algos
	static double Area(const ContourTree& c);
	
	//returns true if c1 and c2 have common area (not a point, but may be an edge)
	static bool DoIntersect(const ContourTree& t1, const Contour& c2);

	//remove points which lie on the same edge
	//removes zero length edges from all contours
	static ContourTree Simplified(const ContourTree& t1);

	//If tree is ok returns true.
	//contacts are allowed
	static bool CheckNoCross(const ContourTree& c);
	//contacts are forbidden
	static bool CheckNoContact(const ContourTree& c);
};

// ContourTree + set of open contours
struct ExtendedTree: public ContourTree {
	ShpVector<Contour> open_contours;

	Contour* get_contour(Point* p) const override;
	Contour* get_contour(Edge* p) const override;
	Contour* get_contour(int i) const override;
	int cont_count() const override { return nodes.size() + open_contours.size(); }
	
	//Methods
	void AddContour(shared_ptr<Contour>& c) override;
	void AddOpenContour(shared_ptr<Contour>& c);
	void RemovePoints(const vector<const Point*>& p) override;
	
	void ReloadEdges() override;

	//Algos
	static ExtendedTree Simplified(const ExtendedTree& t1);
	static ExtendedTree Assemble(const ECollection&);
	static ContourTree ExtractTree(const ExtendedTree& et){ return ContourTree(et); }
};




}
#endif
