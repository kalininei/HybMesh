#include "surface_tree.hpp"

using namespace HM3D;
using namespace HM3D::Surface;

Tree Tree::Assemble(const FaceData& idata){
	// !!! only bounding box check here
	vector<FaceData> data = HM3D::SplitData(idata);
	vector<FaceData*> closed;
	vector<BoundingBox3D> closed_bboxes;
	vector<FaceData*> open;
	//find all closed contours, assemble bounding boxes for them
	for (int i=0; i<data.size(); ++i){
		bool isclosed = IsClosed(data[i]);
		if (isclosed) closed.push_back(&data[i]);
		else open.push_back(&data[i]);
		
		auto ae = AllEdges(data[i]);
		if (isclosed){
			HM3D::VertexData av;
			for (auto e: ae) aa::constant_ids_pvec(e->vertices, 0);
			for (auto e: ae) for (auto v: e->vertices) if (v->id==0){
				av.push_back(v);
				v->id = 1;
			}
			closed_bboxes.push_back(BoundingBox3D(av));
		}
	}
	//add open contours to result
	Tree ret;
	for (auto s: open){ ret.nodes.emplace_back(new Tree::TNode(*s, -1)); }
	if (closed.size() == 0) return ret;
	//add closed contours
	for (auto s: closed) ret.nodes.emplace_back(new Tree::TNode(*s, 0));
	//build tree
	for (int i=1; i<closed.size(); ++i){
		vector<int> inners, outers;
		for (int j=0; j<i; ++j){
			int pos = closed_bboxes[i].relation(closed_bboxes[j]);
			if (pos == 1){
				inners.push_back(j);
			} else if (pos == 2){
				outers.push_back(j);
			} else if (pos == 3){
				//just a sibling, do nothing
			} else throw std::runtime_error(
				"Failed to assemble surface tree for "
				"complicated geometry");
		}
		//look up
		int lowest_outer=-1;
		int klowest = -1;
		for (int k: outers){
			int olv = ret.nodes[open.size()+k]->level;
			if (olv > lowest_outer){
				lowest_outer = olv;
				klowest = k;
			}
		}
		if (klowest >= 0){
			auto out = ret.nodes[open.size() + klowest];
			auto self = ret.nodes[open.size() + i];
			out->children.push_back(self);
			self->parent = out;
			self->level = lowest_outer+1;
		}
		//look down
		for (int k: inners){
			auto in = ret.nodes[open.size() + k];
			auto self = ret.nodes[open.size() + i];
			++(in->level);
			if (in->level == self->level + 1){
				if (!in->parent.expired()){
					auto par = in->parent.lock();
					auto fnd = std::find_if(par->children.begin(), par->children.end(),
							[&in](weak_ptr<Tree::TNode> e){
								return (e.lock() == in);
							});
					assert(fnd != par->children.end());
					par->children.erase(fnd);
				}
				in->parent = self;
				self->children.push_back(in);
			}
		}
	}
	return ret;
}

vector<Tree> Tree::CropLevel01(const Tree& tree){
	vector<Tree> ret;
	for (auto nd: tree.nodes) if (nd->isouter()) {
		ret.push_back(Tree());
		auto& st = ret.back();
		st.nodes.reserve(1 + nd->children.size());
		st.nodes.emplace_back(new Tree::TNode(nd->surface, 0));
		auto& sr = st.nodes[0];
		for (auto ch: nd->children){
			st.nodes.emplace_back(new Tree::TNode(ch.lock()->surface, 1));
			auto& sch = st.nodes.back();
			sch->parent = sr;
			sr->children.push_back(sch);
		}
	}
	return ret;
}
ShpVector<Tree::TNode> Tree::roots() const {
	ShpVector<TNode> ret;
	for (int n=0; n<nodes.size(); ++n){
		if (nodes[n]->level == 0) ret.push_back(nodes[n]);
	}
	return ret;
}

HM3D::FaceData Tree::allfaces() const{
	FaceData ret;
	for (auto& n: nodes){
		ret.insert(ret.end(), n->surface.begin(), n->surface.end());
	}
	return ret;
}

