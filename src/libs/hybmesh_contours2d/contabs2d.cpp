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
vector<VertexEdgeR> Connectivity::VertexEdge(const EdgeData& data, const VertexData& vdata){
	for (auto e: data){
		e->vertices[0]->id = -1;
		e->vertices[1]->id = -1;
	}
	vector<VertexEdgeR> ret(vdata.size());
	for (int i=0; i<vdata.size(); ++i){
		ret[i].v = vdata[i];
		ret[i].v->id = i;
	}
	for (int i=0; i<data.size(); ++i){
		int id1 = data[i]->first()->id;
		int id2 = data[i]->last()->id;
		if (id1>=0) ret[id1].eind.push_back(i);
		if (id2>=0) ret[id2].eind.push_back(i);
	}
	return ret;
}
namespace{
double ameas(double y, double x){
	//this is a faster substitution for atan2 function
	//returns value in [0, 4] ascending by angle.
	double r = (y*y)/(x*x + y*y); //sqr(sin)
	if (x>=0 && y>=0){
		return r;
	} else if (x<0 && y>=0){
		return 2. - r;
	} else if (x<0 && y<0){
		return 2. + r;
	} else {
		return 4. - r;
	}
}

void sort_eind(VertexEdgeR& vr, const EdgeData& data){
	if (vr.size()<3) return;
	//get vectors end points
	vector<Point*> p2(vr.size());
	for (int i=0; i<vr.size(); ++i){
		auto& e = data[vr.eind[i]];
		p2[i] = e->sibling(vr.v.get()).get();
	}
	//calculate angle measures
	vector<double> meas(vr.size());
	for (int i=0; i<vr.size(); ++i){
		meas[i] = ameas(p2[i]->y-vr.v->y, p2[i]->x-vr.v->x);
	}
	//sort measures by value
	vector<double*> pmeas(meas.size());
	for (int i=0; i<meas.size(); ++i) pmeas[i] = &meas[i];
	std::sort(pmeas.begin(), pmeas.end(), [](double* a, double* b){ return *a<*b; });
	//rebuild eind vector
	vector<int> sorted_eind(vr.size());
	for (int i=0; i<vr.size(); ++i){
		int ind = pmeas[i] - &meas[0];
		sorted_eind[i] = vr.eind[ind];
	}
	std::swap(sorted_eind, vr.eind);
}

};

vector<VertexEdgeR> Connectivity::VertexEdgeSorted(const EdgeData& data, const VertexData& vdata){
	vector<VertexEdgeR> ret = VertexEdge(data, vdata);
	for(auto& r: ret) sort_eind(r, data);
	return ret;
}

vector<VertexEdgeR> Connectivity::VertexEdgeSorted(const EdgeData& data){
	vector<VertexEdgeR> ret = VertexEdge(data);
	for(auto& r: ret) sort_eind(r, data);
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

