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

	//Methods
	virtual void AddContour(shared_ptr<Contour>& c);

	//using node->parents calculate children + checks direction
	void UpdateTopology();

	//Algos
	static double Area(const ContourTree& c);
};

// ContourTree + set of open contours
struct ExtendedTree: public ContourTree {
	ShpVector<Contour> open_contours;

	Contour* get_contour(Point* p) const override;
	int cont_count() const override { return nodes.size() + open_contours.size(); }
	
	//Methods
	void AddContour(shared_ptr<Contour>& c) override;

	//Algos
	static ExtendedTree Assemble(const ECollection&);
};




}
#endif
