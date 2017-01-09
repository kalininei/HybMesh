#include "tree.hpp"
#include "cont_assembler.hpp"

using namespace HM2D;
using namespace HM2D::Contour;

double Tree::area() const{
	double ret = 0;
	for (auto& n: nodes) if (n->isbound()){
		double a = fabs(Contour::Area(n->contour));
		if (n->level % 2 == 0) ret += a;
		else ret -= a;
	}
	return ret;
};

Tree Tree::Assemble(const EdgeData& input){
	Tree ret;
	//1) assemble all possible contours
	std::vector<EdgeData> conts = Assembler::AllContours(input);
	//2) add them to ret
	std::for_each(conts.begin(), conts.end(), [&ret](EdgeData& a){ ret.add_contour(a); });
	return ret;
}

EdgeData Tree::alledges() const {
	EdgeData ret;
	for (auto& n: nodes){
		ret.insert(ret.end(), n->contour.begin(), n->contour.end());
	}
	return ret;
};

void Tree::remove_detached(){
	auto s = std::remove_if(nodes.begin(), nodes.end(), [&](shared_ptr<TNode> nd){
			return nd->isdetached(); });
	nodes.resize(s - nodes.begin());
}

ShpVector<Tree::TNode> Tree::detached_contours() const{
	ShpVector<TNode> ret;
	for (auto n: nodes) if (n->isdetached()) ret.push_back(n);
	return ret;
}
ShpVector<Tree::TNode> Tree::bound_contours() const{
	ShpVector<TNode> ret;
	for (auto n: nodes) if (n->isbound()) ret.push_back(n);
	return ret;
}

ShpVector<Contour::Tree::TNode> Contour::Tree::roots() const{
	ShpVector<TNode> ret;
	std::for_each(nodes.begin(), nodes.end(), 
		[&](shared_ptr<TNode> nd){ if (nd->level==0) ret.push_back(nd); }
	);
	return ret;
}

namespace{

//returns whether node was embedded
void add_contour_recursive(shared_ptr<Tree::TNode> node,
		ShpVector<Tree::TNode>& level, shared_ptr<Tree::TNode> parent){
	Point p = InnerPoint(node->contour);
	for (auto& lv: level){
		bool node_within_lv = fabs(Area(lv->contour)) > fabs(Area(node->contour));
		if (WhereIs(lv->contour, p) != INSIDE) node_within_lv = false;
		if (node_within_lv) {
			ShpVector<Tree::TNode> newlevel;
			for (auto w: lv->children) newlevel.push_back(w.lock());
			return add_contour_recursive(node, newlevel, lv);
		}
	}
	node->parent = parent;
	for (auto& lv: level){
		Point p2 = InnerPoint(lv->contour);
		bool lv_within_node = fabs(Area(lv->contour)) < fabs(Area(node->contour));
		if (WhereIs(node->contour, p2) != INSIDE) lv_within_node = false;
		if (lv_within_node) lv->parent = node;
	}
}

}

void Tree::add_contour(const EdgeData& c){
	if (IsOpen(c)){
		return add_detached_contour(c);
	}
	//create node
	auto node = std::make_shared<TNode>(c, 0);
	//call recursive algo
	auto rs = roots();
	add_contour_recursive(node, rs, 0);
	nodes.push_back(node);
	update_topology();
}

void Tree::add_detached_contour(const EdgeData& cont){
	nodes.emplace_back(new TNode(cont, -1));
}

void Contour::Tree::update_topology(){
	//clear children array
	std::for_each(nodes.begin(), nodes.end(),
		[&](shared_ptr<TNode> a){ a->children.clear(); });
	//restore children
	for (auto& c: nodes){
		if (!c->parent.expired()) c->parent.lock()->children.push_back(c);
	}
	//levels
	for (auto& n: nodes) if (IsClosed(n->contour)){
		n->level = 0;
		auto n2 = n;
		while (!n2->parent.expired()){
			++n->level;
			n2 = n2->parent.lock();
		}
	} else { n->level = -1; }
}

Tree Tree::DeepCopy(const Tree& from, int level){
	Tree ret;
	ret.nodes = from.nodes;
	if (level == 0) return ret;
	//copy tree nodes
	for (int i=0; i<ret.nodes.size(); ++i){
		ret.nodes[i] = std::make_shared<TNode>(*ret.nodes[i]);
	}
	std::map<TNode*, int> nind;
	for (int i=0; i<from.nodes.size(); ++i) nind.emplace(from.nodes[i].get(), i);
	for (int i=0; i<ret.nodes.size(); ++i){
		if (!ret.nodes[i]->parent.expired()){
			int ind = nind[ret.nodes[i]->parent.lock().get()];
			ret.nodes[i]->parent = ret.nodes[ind];
		}
		for (auto& w: ret.nodes[i]->children){
			int ind = nind[w.lock().get()];
			w = ret.nodes[ind];
		}
	}
	if (level == 1) return ret;
	//copy contours
	for (int i=0; i<ret.nodes.size(); ++i){
		EdgeData newed;
		HM2D::DeepCopy(ret.nodes[i]->contour, newed, level-2);
		std::swap(ret.nodes[i]->contour, newed);
	}
	return ret;
}

shared_ptr<Tree::TNode> Tree::find_node(const Point* p) const{
	shared_ptr<Tree::TNode> ret;
	for (auto n: nodes){
		if (Contains(n->contour, p)){
			ret = n;
			break;
		}
	}
	return ret;
}
shared_ptr<Tree::TNode> Tree::find_node(const Edge* p) const{
	shared_ptr<Tree::TNode> ret;
	for (auto n: nodes){
		if (Contains(n->contour, p)){
			ret = n;
			break;
		}
	}
	return ret;
}

int Tree::whereis(const Point& p) const{
	//using WhereIs routine for single contour
	//since algorithms are the same:
	//calculation of crosses between [p, far point]
	EdgeData cedges;
	for (auto n: nodes) if (n->isbound()){
		cedges.insert(cedges.end(), n->contour.begin(),
				n->contour.end());
	}
	return Contour::WhereIs(cedges, p);
}
