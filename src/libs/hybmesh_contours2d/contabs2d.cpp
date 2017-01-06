#include "contabs2d.hpp"
using namespace HM2D;
using namespace HM2D::Connectivity;

vector<VertexEdgeR> Connectivity::VertexEdge(const EdgeData& data){
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

vector<vector<int>> Connectivity::EdgeEdge(const EdgeData& data){
	auto ve = VertexEdge(data);
	vector<vector<int>> ret(data.size());
	for (auto& c: ve){
		for (int i1=0; i1<c.size(); ++i1)
		for (int i2=i1+1; i2<c.size(); ++i2){
			ret[c.eind[i1]].push_back(c.eind[i2]);
			ret[c.eind[i2]].push_back(c.eind[i1]);
		}
	}

	return ret;
}
