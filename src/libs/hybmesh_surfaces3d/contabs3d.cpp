#include "contabs3d.hpp"

namespace ct = HM3D::Connectivity;
using namespace ct;

vector<VertexEdgeR> ct::VertexEdge(const EdgeData& data){
	auto ev = AllVertices(data);
	vector<VertexEdgeR> ret(ev.size());
	for (int i=0; i<ev.size(); ++i){
		ret[i].v = ev[i];
		ret[i].v->id = i;
	}
	for (int i=0; i<data.size(); ++i){
		ret[data[i]->first()->id].eind.push_back(i);
		ret[data[i]->last()->id].eind.push_back(i);
	}
	return ret;
}

vector<EdgeFaceR> ct::EdgeFace(const FaceData& data){
	auto ae = AllEdges(data);
	vector<EdgeFaceR> ret(ae.size());
	for (int i=0; i<ae.size(); ++i){
		ret[i].e = ae[i];
		ret[i].e->id = i;
	}
	for (int i=0; i<data.size(); ++i)
	for (auto& e: data[i]->edges){
		ret[e->id].find.push_back(i);
	}
	return ret;
}

vector<EdgeFaceExtendedR> ct::EdgeFaceExtended(const FaceData& data){
	auto ae = AllEdges(data);
	vector<EdgeFaceExtendedR> ret(ae.size());
	for (int i=0; i<ae.size(); ++i){
		ret[i].e = ae[i];
		ret[i].e->id = i;
	}
	for (int i=0; i<data.size(); ++i)
	for (int j=0; j<data[i]->edges.size(); ++j){
		auto& r = ret[data[i]->edges[j]->id];
		r.find.push_back(i);
		r.locind.push_back(j);
		r.posdir.push_back(data[i]->is_positive_edge(j));
	}
	return ret;
}

vector<vector<int>> ct::FaceFace(const ShpVector<Face>& data){
	return FaceFace(EdgeFace(data), data.size());
}

vector<vector<int>> ct::FaceFace(const vector<EdgeFaceR>& edge_face, int nfaces){
	vector<vector<int>> ret(nfaces);
	for (auto& it: edge_face){
		for (int i=0; i<it.size(); ++i){
			int f1 = it.find[i];
			for (int j=i+1; j<it.size(); ++j){
				int f2 = it.find[j];
				ret[f1].push_back(f2);
				ret[f2].push_back(f1);
			}
		}
	}
	return ret;
}

vector<vector<int>> ct::CellCell(const CellData& data){
	vector<vector<int>> ret(data.size());
	auto af = AllFaces(data);
	aa::enumerate_ids_pvec(data);
	for (auto f: af) if (!f->is_boundary()){
		auto c1 = f->left.lock()->id;
		auto c2 = f->right.lock()->id;
		ret[c1].push_back(c2);
		ret[c2].push_back(c1);
	}
	return ret;
}
