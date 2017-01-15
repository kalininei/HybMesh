#include "contabs2d.hpp"
#include "contour.hpp"
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

vector<vector<int>> Connectivity::CellCell(const CellData& data){
	vector<vector<int>> ret(data.size());
	auto ae = AllEdges(data);
	aa::enumerate_ids_pvec(data);
	for (auto e: ae) if (!e->is_boundary()){
		auto c1 = e->left.lock()->id;
		auto c2 = e->right.lock()->id;
		ret[c1].push_back(c2);
		ret[c2].push_back(c1);
	}

	return ret;
}

vector<Connectivity::VertexCellR> Connectivity::VertexCell(const CellData& data){
	auto av = AllVertices(data);
	vector<Connectivity::VertexCellR> ret(av.size());
	for (int i=0; i<av.size(); ++i) ret[i].v = av[i];
	aa::enumerate_ids_pvec(av);
	for (int i=0; i<data.size(); ++i){
		auto op = Contour::OrderedPoints(data[i]->edges);
		for (int j=0; j<op.size()-1; ++j){
			ret[op[j]->id].cind.push_back(i);
		}
	}
	return ret;
}

