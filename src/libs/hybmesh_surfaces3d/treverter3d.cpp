#include <stack>
#include "treverter3d.hpp"

using namespace HM3D;
using namespace HM3D::Surface;
using namespace HM3D::Surface::R;

Revert::Revert(const FaceData& srf){
	permanent = false;
	obj = const_cast<FaceData*>(&srf);
	need_revert = vector<bool>(srf.size(), false);
	auto ae = AllEdges(srf);
	aa::enumerate_ids_pvec(ae);
	aa::enumerate_ids_pvec(srf);
	vector<int> adjfaces1(ae.size(), -1), adjfaces2(ae.size(), -1);
	vector<int> adjlocind1(ae.size(), -1), adjlocind2(ae.size(), -1);
	for (size_t i=0; i<srf.size(); ++i)
	for (size_t j=0; j<srf[i]->edges.size(); ++j){
		int eid = srf[i]->edges[j]->id;
		if (adjfaces1[eid] == -1){
			adjfaces1[eid] = i;
			adjlocind1[eid] = j;
		} else if (adjfaces2[eid] == -1){
			adjfaces2[eid] = i;
			adjlocind2[eid] = j;
		} else throw std::runtime_error("edge has more than two adjacent faces");
	}
	vector<bool> processed_faces(srf.size(), false);
	//edge<0> in face<1> should be positive/negative<2>
	std::stack<std::tuple<int, int, bool>> face_list;

	while (1){
		auto fnd = std::find(processed_faces.begin(), processed_faces.end(), false);
		if (fnd == processed_faces.end()) break;
		face_list.push(std::make_tuple(fnd - processed_faces.begin(), 0, true));
		*fnd = true;
		while (face_list.size() > 0){
			auto& finfo = face_list.top();
			int& iedge = std::get<0>(finfo);
			bool& should_be_positive = std::get<2>(finfo);
			auto fc = srf[std::get<1>(finfo)];
			bool really_positive = fc->is_positive_edge(iedge);
			if (should_be_positive != really_positive){
				need_revert[fc->id] = true;
			}
			face_list.pop();
			for (int j=0; j<fc->edges.size(); ++j){
				int ied = fc->edges[j]->id;
				int afc=-1, locind;
				if (adjfaces1[ied] == fc->id){
					afc = adjfaces2[ied];
					locind = adjlocind2[ied];
				} else {
					assert(adjfaces2[ied] == fc->id);
					afc = adjfaces1[ied];
					locind = adjlocind1[ied];
				}
				if (afc >= 0 && processed_faces[afc] == false){
					bool pos = fc->is_positive_edge(j);
					if (need_revert[fc->id]) pos = !pos;
					processed_faces[afc] = true;
					face_list.push(std::make_tuple(locind, afc, !pos));
				}
			}
		}
	}
	//revert
	for (size_t i=0; i<obj->size(); ++i)
	if (need_revert[i]){
		(*obj)[i]->reverse();
	}
}
Revert::~Revert(){
	if (!permanent){
		for (size_t i=0; i<obj->size(); ++i)
		if (need_revert[i]){
			(*obj)[i]->reverse();
		}
	}
}
void Revert::reverse_direction(){
	for (size_t i=0; i<obj->size(); ++i){
		need_revert[i] = !need_revert[i];
		(*obj)[i]->reverse();
	}
}

RevertTree::RevertTree(const Tree& srf){
	obj = const_cast<Tree*>(&srf);
	for (auto& n: obj->nodes){
		if (n->isdetached()) openrevs.emplace_back(new Revert(n->surface));
		else closedrevs.emplace_back(new Revert(n->surface));
	}
	int k=0;
	for (auto& n: obj->nodes) if (n->isbound()){
		auto& r = closedrevs[k++];
		double volume = Surface::Volume(n->surface);
		if (volume < 0 && n->level % 2 == 0) r->reverse_direction();
		else if (volume > 0 && n->level % 2 == 1) r->reverse_direction();
	}
}
RevertTree::~RevertTree(){}

RevertGridSurface::RevertGridSurface(const FaceData& srf, bool cells_left){
	permanent = false;
	obj = const_cast<FaceData*>(&srf);
	need_revert.resize(srf.size(), false);
	if (cells_left) for (int i=0; i<srf.size(); ++i){
		auto& f = srf[i];
		if (!f->has_left_cell() && f->has_right_cell()){
			need_revert[i] = true;
			f->reverse();
		}
	} else for (int i=0; i<srf.size(); ++i){
		auto& f = srf[i];
		if (!f->has_right_cell() && f->has_left_cell()){
			need_revert[i] = true;
			f->reverse();
		}
	}
}
RevertGridSurface::~RevertGridSurface(){
	if (!permanent){
		for (int i=0; i<obj->size(); ++i) if (need_revert[i]){
			(*obj)[i]->reverse();
		}
	}
}
