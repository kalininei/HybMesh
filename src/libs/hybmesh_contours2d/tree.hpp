#ifndef HMCONT2D_TREE_HPP
#define HMCONT2D_TREE_HPP

#include "collections.hpp"
#include "contour.hpp"

namespace HMCont2D{
struct ContourTree;
struct ExtendedTree;

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
	virtual void AddContour(const Contour& c);
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
	//build an extended tree with 0 open contours
	static ExtendedTree AsExtended(const ContourTree& et);
	//vtk save
	static void SaveVtk(const ContourTree& ct, const char* fn);

	//If tree is ok returns true.
	//contacts are allowed
	static bool CheckNoCross(const ContourTree& c);
	//contacts are forbidden
	static bool CheckNoContact(const ContourTree& c);

	//DeepCopy overriden from ECollection
	template<class TTarget, class = Tpp::IsBase<ContourTree, TTarget>>
	static TDeepCopyResult DeepCopy(const ContourTree& from, TTarget& to, const ShpGen& gen) {
		TDeepCopyResult res;
		//copy nodes
		for (auto c: from.nodes){
			//to.nodes.push_back(TreeNode());
			aa::add_shared(to.nodes, TreeNode());
			TDeepCopyResult dcr = Contour::DeepCopy(*c, *to.nodes.back(), gen);
			res.oldnew.insert(dcr.oldnew.begin(), dcr.oldnew.end());
			res.newold.insert(dcr.newold.begin(), dcr.newold.end());
		}
		//fill parent and childs
		for (int i=0; i<from.nodes.size(); ++i){
			auto par = from.nodes[i]->parent;
			if (par == NULL) to.nodes[i]->parent = 0;
			else{
				for (int j=0; j<from.nodes.size(); ++j){
					if (from.nodes[j].get() == par){
						to.nodes[i]->parent = to.nodes[j].get();
						break;
					}
				}
			}
		}
		to.UpdateTopology();
		//fill data
		to.ReloadEdges();
		return res;
	};
};

// ContourTree + set of open contours
struct ExtendedTree: public ContourTree {
	ShpVector<Contour> open_contours;

	Contour* get_contour(Point* p) const override;
	Contour* get_contour(Edge* p) const override;
	Contour* get_contour(int i) const override;
	vector<Contour*> all_contours() const;
	int cont_count() const override { return nodes.size() + open_contours.size(); }
	
	//Methods
	void AddContour(shared_ptr<Contour>& c) override;
	void AddContour(const Contour& c) override;
	void AddOpenContour(shared_ptr<Contour>& c);
	void RemovePoints(const vector<const Point*>& p) override;
	void ReloadEdges() override;

	//Algos
	static ContourTree ExtractTree(const ExtendedTree& et){ return ContourTree(et); }
	//vtk save
	static void SaveVtk(const ExtendedTree& ct, const char* fn);
	
	//DeepCopy overriden from ECollection
	template<class TTarget, class = Tpp::IsBase<ExtendedTree, TTarget>>
	static TDeepCopyResult DeepCopy(const ExtendedTree& from, TTarget& to, const ShpGen& gen) {
		TDeepCopyResult res = ContourTree::DeepCopy(from, to, gen);
		//copy open contours
		for (auto c: from.open_contours){
			aa::add_shared(to.open_contours, Contour());
			TDeepCopyResult dcr = Contour::DeepCopy(*c, *to.open_contours.back(), gen);
			res.oldnew.insert(dcr.oldnew.begin(), dcr.oldnew.end());
			res.newold.insert(dcr.newold.begin(), dcr.newold.end());
		}
		//fill data
		to.ReloadEdges();
		return res;
	};
};




}
#endif
