#include "tree.hpp"
#include "contclipping.hpp"
#include <fstream>

using namespace HMCont2D;

int ContourTree::TreeNode::level() const{
	if (parent == NULL) return 0;
	else return parent->level() + 1;
}

vector<ContourTree::TreeNode*> ContourTree::roots() const{
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

ExtendedTree ContourTree::AsExtended(const ContourTree& et){
	ExtendedTree ret;
	ret.data = et.data;
	ret.nodes = et.nodes;
	return ret;
}

Contour* ContourTree::get_contour(Point* p) const{
	Contour* ret = 0;
	auto fnd = std::find_if(nodes.begin(), nodes.end(),
		[&p](shared_ptr<Contour> a){ return a->contains_point(p); }
	);
	if (fnd != nodes.end()) return fnd->get();
	else return 0;
}

Contour* ContourTree::get_contour(Edge* p) const{
	Contour* ret = 0;
	auto fnd = std::find_if(nodes.begin(), nodes.end(),
		[&p](shared_ptr<Contour> a){ return a->contains(p); }
	);
	if (fnd != nodes.end()) return fnd->get();
	else return 0;
}

Contour* ExtendedTree::get_contour(Edge* e) const{
	Contour* ret = ContourTree::get_contour(e);
	if (ret != 0) return ret;
	for (auto& oc: open_contours) if (oc->contains(e)) return oc.get();
	return 0;
}

vector<Contour*> ExtendedTree::all_contours() const{
	vector<Contour*> ret;
	for (auto n: nodes) ret.push_back(n.get());
	for (auto n: open_contours) ret.push_back(n.get());
	return ret;
}

namespace{

//returns whether node was embedded
void add_contour_recursive(ContourTree::TreeNode* node,
		vector<ContourTree::TreeNode*>& level, ContourTree::TreeNode* parent){
	Point p = node->InnerPoint();
	for (auto& lv: level){
		bool node_within_lv = lv->IsWithin(p) &&
			fabs(HMCont2D::Contour::Area(*lv)) > fabs(HMCont2D::Contour::Area(*node));
		if (node_within_lv) return add_contour_recursive(node, lv->children, lv);
	}
	node->parent = parent;
	for (auto& lv: level){
		Point p2 = lv->InnerPoint();
		bool lv_within_node = node->IsWithin(p2) &&
			fabs(HMCont2D::Contour::Area(*lv)) < fabs(HMCont2D::Contour::Area(*node));
		if (lv_within_node) lv->parent = node;
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
	add_contour_recursive(node.get(), rs, 0);
	nodes.push_back(node);
	UpdateTopology();
}

void ContourTree::AddContour(const Contour& c){
	shared_ptr<Contour> c2(new Contour(c));
	AddContour(c2);
}

void ContourTree::RemovePoints(const vector<const Point*>& p){
	for (auto& n: nodes) n->RemovePoints(p);
	ReloadEdges();
}

void ContourTree::ReloadEdges(){
	data.clear();
	for (auto& n: nodes) data.insert(data.end(), n->data.begin(), n->data.end());
}
void ContourTree::Reallocate(){
	for (auto& n: nodes) n->Reallocate();
	ReloadEdges();
}

void ExtendedTree::RemovePoints(const vector<const Point*>& p){
	for (auto& n: nodes) n->RemovePoints(p);
	for (auto& n: open_contours) n->RemovePoints(p);
	ReloadEdges();
}

void ExtendedTree::ReloadEdges(){
	data.clear();
	for (auto& n: nodes) data.insert(data.end(), n->data.begin(), n->data.end());
	for (auto& n: open_contours) data.insert(data.end(), n->data.begin(), n->data.end());
}

void ExtendedTree::Reallocate(){
	for (auto& n: nodes) n->Reallocate();
	for (auto& n: open_contours) n->Reallocate();
	ReloadEdges();
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

bool ContourTree::IsWithin(const Point& p) const{
	for (auto nd: nodes) if (HMCont2D::Contour::Area(*nd) > 0){
		if (nd->IsWithin(p)){
			for (auto ch: nd->children){
				if (ch->IsWithin(p)) goto NEXTNODE;
			}
			return true;
		}
	NEXTNODE:;
	}
	return false;
}

bool ContourTree::IsWithout(const Point& p) const{
	for (auto nd: nodes) if (HMCont2D::Contour::Area(*nd) > 0){
		if (!nd->IsWithout(p)){
			for (auto ch: nd->children){
				if (ch->IsWithin(p)) goto NEXTNODE;
			}
			return false;
		}
	NEXTNODE:;
	}
	return true;
}


void ExtendedTree::AddContour(shared_ptr<Contour>& c){
	if (c->is_closed()) ContourTree::AddContour(c);
	else AddOpenContour(c);
}

void ExtendedTree::AddContour(const Contour& c){
	if (c.is_closed()) ContourTree::AddContour(c);
	else {
		auto sp = std::make_shared<Contour>(c);
		AddOpenContour(sp);
	}
}

void ExtendedTree::AddOpenContour(shared_ptr<Contour>& c){
	open_contours.push_back(c);
	Unite(*c);
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


namespace {

bool check_core(const ContourTree& tree, double v0, double v1){
	//all contours are closed
	for (auto cont: tree.nodes){
		if (!cont->is_closed()) return false;
	}
	//each edge is owned by one contour
	std::set<Edge*> eset;
	for (auto cont: tree.nodes){
		for (auto e: cont->data){
			if (eset.find(e.get()) != eset.end()) return false;
			eset.insert(e.get());
		}
	}
	//no edge crosses
	double ksieta[2];
	for (int i=0; i<tree.size(); ++i){
		for (int j=i+1; j<tree.size(); ++j){
			Edge* e1 = tree.edge(i);
			Edge* e2 = tree.edge(j);
			if (Edge::AreConnected(*e1, *e2)) continue;
			SectCross(*e1->pstart, *e1->pend, *e2->pstart, *e2->pend, ksieta);
			if (ksieta[0]>v0 && ksieta[0]<v1 && ksieta[1]>v0 && ksieta[1]<v1)
				return false;
		}
	}
	return true;
}

}

bool ContourTree::CheckNoCross(const ContourTree& tree){
	return check_core(tree, +geps, 1-geps);
}

bool ContourTree::CheckNoContact(const ContourTree& tree){
	return check_core(tree, -geps, 1+geps);
}


Contour* ContourTree::get_contour(int i) const{
	return nodes[i].get();
}

Contour* ExtendedTree::get_contour(int i) const{
	if (i < ContourTree::cont_count()) return ContourTree::get_contour(i);
	else{
		i = i - ContourTree::cont_count();
		return open_contours[i].get();
	}
}

void ContourTree::SaveVtk(const ContourTree& ct, const char* fn){
	auto et = AsExtended(ct);
	ExtendedTree::SaveVtk(et, fn);
}

void ExtendedTree::SaveVtk(const ExtendedTree& ct, const char* fn){
	//save basic collection
	ECollection::SaveVtk(ct, fn);
	//set contours ids, local edge indicies
	vector<int> cont_ids(ct.data.size(), -1);
	vector<int> loc_ids(ct.data.size(), -1);
	//closed contours
	int ind = 0;
	for (auto cont: ct.nodes){
		auto locind = 0;
		for (auto e: *cont){
			int gind = ct.get_index(e.get());
			cont_ids[gind] = ind;
			loc_ids[gind] = locind++;
		}
		++ind;
	}
	//open contours
	for (auto cont: ct.open_contours){
		auto locind = 0;
		for (auto e: *cont){
			int gind = ct.get_index(e.get());
			cont_ids[gind] = ind;
			loc_ids[gind] = locind++;
		}
		++ind;
	}
	
	//write to file
	std::ofstream of(fn, std::ios::app);
	of<<"CELL_DATA "<<cont_ids.size()<<std::endl;
	of<<"SCALARS contour_id int 1"<<std::endl;
	of<<"LOOKUP_TABLE default"<<std::endl;
	for (auto i: cont_ids) of<<i<<std::endl;
	of<<"SCALARS local_edge_index int 1"<<std::endl;
	of<<"LOOKUP_TABLE default"<<std::endl;
	for (auto i: loc_ids) of<<i<<std::endl;

}
