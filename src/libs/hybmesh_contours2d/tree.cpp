#include "tree.hpp"

using namespace HMCont2D;

vector<ContourTree::TreeNode*> ContourTree::roots(){
	vector<TreeNode*> ret;
	std::for_each(nodes.begin(), nodes.end(), 
		[&](shared_ptr<TreeNode> nd){ if (nd->parent == 0) ret.push_back(nd.get()); }
	);
	return ret;
}

double ContourTree::Area(const ContourTree& c){
	return std::accumulate(c.nodes.begin(), c.nodes.end(), 0.0,
			[](double s, shared_ptr<Contour> cc){ return s + Contour::Area(*cc); });
}

Contour* ContourTree::get_contour(Point* p) const{
	Contour* ret = 0;
	auto fnd = std::find_if(nodes.begin(), nodes.end(),
		[&p](shared_ptr<Contour> a){ return a->contains_point(p); }
	);
	if (fnd != nodes.end()) return fnd->get();
	else return 0;
}

namespace{

//returns whether node was embedded
void AddContourRecursive(ContourTree::TreeNode* node,
		vector<ContourTree::TreeNode*>& level, ContourTree::TreeNode* parent){
	Point* p = node->first();
	for (auto& lv: level){
		if (lv->IsWithin(*p)){
			return AddContourRecursive(node, lv->children, lv);
		}
	}
	node->parent = parent;
	for (auto& lv: level){
		Point* p2 = lv->first();
		if (node->IsWithin(*p2)){
			lv->parent = node;
		}
	}
}

}

void ContourTree::AddContour(shared_ptr<Contour>& c){
	assert(c->is_closed());
	//create node
	shared_ptr<TreeNode> node(new TreeNode);
	node->Unite(*c);
	c = node;
	Unite(*node);
	//call recursive algo
	auto rs = roots();
	AddContourRecursive(node.get(), rs, 0);
	nodes.push_back(node);
	UpdateTopology();
}

void ContourTree::UpdateTopology(){
	//clear children array
	std::for_each(nodes.begin(), nodes.end(),
		[&](shared_ptr<TreeNode> a){ a->children.clear(); });
	//restore children
	for (auto& c: nodes){
		if (c->parent != 0) c->parent->children.push_back(c.get());
	}
	
	std::function<void(TreeNode*, bool)>
	direct = [&direct](TreeNode* nd, bool dir){
		nd->ForceDirection(dir);
		for (auto& c: nd->children) direct(c, !dir);
	};

	//all roots to positive direction
	for (auto& c: nodes){
		if (c->parent == 0) direct(c.get(), true);
	}
}


void ExtendedTree::AddContour(shared_ptr<Contour>& c){
	if (c->is_closed()){
		ContourTree::AddContour(c);
	} else {
		open_contours.push_back(c);
		Unite(*c);
	}
}

ExtendedTree ExtendedTree::Assemble(const ECollection& col){
	ExtendedTree ret;
	auto ap = col.all_points();
	std::set<Point*> unusedpnts(ap.begin(), ap.end());
	//1) assemble all possible contours
	ShpVector<Contour> conts;
	while (unusedpnts.size() > 0){
		Point* p1 = *unusedpnts.begin();
		conts.push_back(
			std::make_shared<Contour>(Contour::Assemble(col, p1))
		);
		for (Point* p: conts.back()->all_points()){
			unusedpnts.erase(p);
		}
	}
	//2) add'em
	std::for_each(conts.begin(), conts.end(), [&ret](shared_ptr<Contour> a){ ret.AddContour(a); });
	return ret;
}

Contour* ExtendedTree::get_contour(Point* p) const{
	Contour* ret = ContourTree::get_contour(p);
	if (ret == 0){
		auto fnd = std::find_if(open_contours.begin(), open_contours.end(),
			[&p](shared_ptr<Contour> a){ return a->contains_point(p); }
		);
		if (fnd != open_contours.end()) return fnd->get();
	}
	return ret;
}
