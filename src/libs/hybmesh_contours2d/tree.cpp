#include "tree.hpp"
#include "cont_assembler.hpp"
#include "treverter2d.hpp"
#include "contour.hpp"
#include "finder2d.hpp"

using namespace HM2D;
using namespace HM2D::Contour;

double Tree::area() const{
	double ret = 0;
	for (auto& n: nodes) if (n->isbound()){
		double a = fabs(Contour::Area(n->contour));
		if (n->isouter()) ret += a;
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

EdgeData Tree::alledges_bound() const {
	EdgeData ret;
	for (auto& n: bound_contours()){
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
		[&](shared_ptr<TNode> nd){ if (nd->isroot()) ret.push_back(nd); }
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
		if (Contour::Finder::WhereIs(lv->contour, p) != INSIDE) node_within_lv = false;
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
		if (Contour::Finder::WhereIs(node->contour, p2) != INSIDE) lv_within_node = false;
		if (lv_within_node) lv->parent = node;
	}
}

}

void Tree::add_contour(const EdgeData& c){
	EdgeData c2(c);
	add_contour(std::move(c2));
}

void Tree::add_detached_contour(const EdgeData& c){
	EdgeData c2(c);
	add_detached_contour(std::move(c2));
}
void Tree::add_contour(EdgeData&& c){
	if (IsOpen(c)){
		return add_detached_contour(std::move(c));
	}
	//create node
	auto node = std::make_shared<TNode>(std::move(c), 0);
	//call recursive algo
	auto rs = roots();
	add_contour_recursive(node, rs, 0);
	nodes.push_back(node);
	update_topology();
}
void Tree::add_detached_contour(EdgeData&& c){
	nodes.emplace_back(new TNode(std::move(c), -1));
}
void Tree::remove_contour(TNode* n){
	int i = aa::shpvec_ifind(nodes, n);
	assert(i<nodes.size());
	for (auto ch: n->children){
		ch.lock()->parent = n->parent;
	}
	nodes.erase(nodes.begin()+i);
	update_topology();
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
	for (auto& n: nodes) if (n->isbound()){
		n->level = 0;
		auto n2 = n;
		while (!n2->parent.expired()){
			++n->level;
			n2 = n2->parent.lock();
		}
	}
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
		if (HM2D::Finder::Contains(n->contour, p)){
			ret = n;
			break;
		}
	}
	return ret;
}
shared_ptr<Tree::TNode> Tree::find_node(const Edge* p) const{
	shared_ptr<Tree::TNode> ret;
	for (auto n: nodes){
		if (HM2D::Finder::Contains(n->contour, p)){
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
	return Contour::Finder::WhereIs(cedges, p);
}

Tree Tree::GridBoundary(const GridData& data){
	vector<Tree> vt = GridBoundary01(data);
	if (vt.size() == 0) return Tree();
	if (vt.size() == 1) return vt[0];

	Tree ret;
	if (std::all_of(vt.begin(), vt.end(),
			[](const Tree& t){ return t.nodes.size() == 1; })){
		for (int i=0; i<vt.size(); ++i){
			ret.nodes.push_back(vt[i].nodes[0]);
		}
		return ret;
	}

	//add all roots
	for (int i=0; i<vt.size(); ++i){
		ret.nodes.push_back(vt[i].nodes[0]);
	}
	//calculate roots relation
	auto _where_is = [](const TNode& tar, const TNode& src)->int{
		for (auto& e: src.contour){
			int w = HM2D::Contour::Finder::WhereIs(tar.contour, e->center());
			if (w != BOUND) return w;
		}
		return BOUND;
	};
	vector<vector<int>> outsides(ret.nodes.size()), insides(ret.nodes.size());
	for (int i=0; i<ret.nodes.size(); ++i){
		for (int j=i+1; j<ret.nodes.size(); ++j){
			int w = _where_is(*ret.nodes[i], *ret.nodes[j]);
			if (w == INSIDE) {
				outsides[j].push_back(i);
				insides[i].push_back(j);
			} else if (w == OUTSIDE){
				outsides[i].push_back(j);
				insides[j].push_back(i);
			} else {}
		}
	}
	//calculate levels
	for (int i=0; i<ret.nodes.size(); ++i){
		for (auto j: insides[i]) ret.nodes[j]->level+=2;
	}
	//calculate root level parents
	for (int i=0; i<ret.nodes.size(); ++i) if (ret.nodes[i]->level > 0){
		int highest_level = -1;
		for (auto j: outsides[i]){
			if (ret.nodes[j]->level > highest_level){
				ret.nodes[i]->parent = ret.nodes[j];
				highest_level = ret.nodes[j]->level;
			}
		}
	}
	//parents to children nodes
	for (int i=0; i<ret.nodes.size(); ++i) if (ret.nodes[i]->level>0){
		while (!ret.nodes[i]->parent.expired()){
			bool found = false;

			if (ret.nodes[i]->parent.lock()->level != 1)
			for (auto c: ret.nodes[i]->parent.lock()->children){
				if (_where_is(*c.lock(), *ret.nodes[i]) == INSIDE){
					ret.nodes[i]->parent = c;
					found = true;
					break;
				}
			}

			if (found) break;
			else ret.nodes[i]->parent = ret.nodes[i]->parent.lock()->parent;
		}
		assert(ret.nodes[i]->parent.expired() ||
		       ret.nodes[i]->parent.lock()->level == 1);
	}

	//add children and update topology
	for (int i=0; i<vt.size(); ++i){
		ret.nodes.insert(ret.nodes.end(), vt[i].nodes.begin()+1, vt[i].nodes.end());
	}
	ret.update_topology();
	return ret;
}

vector<Tree> Tree::CropLevel01(const Tree& data){
	vector<Tree> ret;
	for (auto& n: data.nodes) if (n->isouter()){
		ret.emplace_back();
		auto& t = ret.back();
		t.nodes.push_back(std::make_shared<TNode>(n->contour, 0));
		for (auto& ch: n->children){
			t.nodes.push_back(std::make_shared<TNode>(
				ch.lock()->contour, 1));
			t.nodes[0]->children.push_back(t.nodes.back());
			t.nodes.back()->parent = t.nodes[0];
		}
	}
	return ret;
}

vector<Tree> Tree::GridBoundary01(const GridData& data){
	vector<EdgeData> conts = HM2D::Contour::Assembler::GridBoundary(data);
	vector<EdgeData*> inners, outers;
	for (auto& c: conts){
		Contour::R::Clockwise rev(c, false);
		if (c[0]->has_left_cell()) outers.push_back(&c);
		else inners.push_back(&c);
	}

	vector<Tree> ret;
	for (int i=0; i<outers.size(); ++i){
		ret.emplace_back();
		ret.back().nodes.push_back(std::make_shared<TNode>(
			std::move(*outers[i]), 0));
	}

	for (auto in: inners){
		std::vector<int> outsiders;
		Point* p = (*in)[0]->vertices[0].get();
		for (int i=0; i<ret.size(); ++i){
			if (Contour::Finder::WhereIs(ret[i].nodes[0]->contour, *p)
					==INSIDE){
				outsiders.push_back(i);
			}
		}
		assert(outsiders.size()>0);
		for (int i=1; i<outsiders.size(); ++i){
			double a1 = ret[outsiders[0]].area();
			double a2 = ret[outsiders[i]].area();
			if (a2 < a1) outsiders[0] = outsiders[i];
		}
		auto& r = ret[outsiders[0]];
		r.nodes.push_back(std::make_shared<TNode>(
			std::move(*in), 1));
		r.nodes[0]->children.push_back(r.nodes.back());
		r.nodes.back()->parent = r.nodes[0];
	}

	return ret;
}
