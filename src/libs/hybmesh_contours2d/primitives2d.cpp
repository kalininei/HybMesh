#include "primitives2d.hpp"
#include "contabs2d.hpp"
#include "hmgraph.hpp"
using namespace HM2D;

void GridData::enumerate_all() const{
	aa::enumerate_ids_pvec(vvert);
	aa::enumerate_ids_pvec(vedges);
	aa::enumerate_ids_pvec(vcells);
}
// ================= Extract primitives
VertexData HM2D::AllVertices(const EdgeData& from){
	for (auto& e: from){
		e->vertices[0]->id = 0;
		e->vertices[1]->id = 0;
	}
	VertexData ret;
	for (auto& e: from){
		if (e->vertices[0]->id == 0){
			e->vertices[0]->id = 1;
			ret.push_back(e->vertices[0]);
		}
		if (e->vertices[1]->id == 0){
			e->vertices[1]->id = 1;
			ret.push_back(e->vertices[1]);
		}
	}
	return ret;
}
VertexData HM2D::AllVertices(const CellData& from){
	return AllVertices(AllEdges(from));
}
EdgeData HM2D::AllEdges(const CellData& from){
	for (auto& c: from)
	for (auto& e: c->edges) e->id = 0;

	EdgeData ret;
	for (auto& c: from)
	for (auto& e: c->edges) if (e->id == 0){
		e->id = 1;
		ret.push_back(e);
	}
	return ret;
}
std::tuple<VertexData> HM2D::AllPrimitives(const EdgeData& from){
	return std::make_tuple(AllVertices(from));
}
std::tuple<VertexData, EdgeData> AllPrimitives(const CellData& from){
	std::tuple<VertexData, EdgeData> ret;
	for (auto& c: from)
	for (auto& e: c->edges){
		e->id = 0;
		e->first()->id = 0;
		e->last()->id = 0;
	}
	for (auto& c: from)
	for (auto& e: c->edges) if (e->id == 0){
		e->id = 1;
		std::get<1>(ret).push_back(e);
		if (e->first()->id == 0){
			e->first()->id = 1;
			std::get<0>(ret).push_back(e->first());
		}
		if (e->last()->id == 0){
			e->last()->id = 1;
			std::get<0>(ret).push_back(e->last());
		}
	}

	return ret;
}


// ================= DeepCopy
namespace {
template<class T>
void shared_vec_deepcopy(const vector<shared_ptr<T>>& from, vector<shared_ptr<T>>& to){
	to.resize(to.size() + from.size());
	auto it1 = to.end() - from.size();
	auto it2 = from.begin();
	while (it1 != to.end()){
		(it1++)->reset(new T(**(it2++)));
	}
}
}
void HM2D::DeepCopy(const VertexData& from, VertexData& to){
	shared_vec_deepcopy(from, to);
}
void HM2D::DeepCopy(const EdgeData& from, EdgeData& to, int level){
	shared_vec_deepcopy(from, to);
	if (level > 0){
		VertexData vorig = AllVertices(from);
		VertexData vnew;
		DeepCopy(vorig, vnew);
		aa::enumerate_ids_pvec(vorig);
		aa::enumerate_ids_pvec(vnew);
		auto eit = to.end() - from.size();
		while (eit != to.end()){
			for (auto& v: (*eit++)->vertices){
				v = vnew[v->id];
			}
		}
	}
}
void HM2D::DeepCopy(const CellData& from, CellData& to, int level){
	shared_vec_deepcopy(from, to);
	if (level > 0){
		EdgeData forig = AllEdges(from);
		EdgeData fnew;
		DeepCopy(forig, fnew, level - 1);
		aa::enumerate_ids_pvec(forig);
		aa::enumerate_ids_pvec(fnew);
		auto cit = to.end() - from.size();
		int k = 0;
		while (cit != to.end()){
			for (auto& e: (*cit)->edges){
				e = fnew[e->id];
				if (forig[e->id]->left.lock() == from[k]) e->left = *cit;
				else if (forig[e->id]->right.lock() == from[k]) e->right = *cit;
			}
			++k;
			++cit;
		}
	}
}
void HM2D::DeepCopy(const GridData& from, GridData& to, int level){
	to.clear();
	if (level>=2) DeepCopy(from.vvert, to.vvert);
	else to.vvert = from.vvert;
	if (level>=1) DeepCopy(from.vedges, to.vedges, 0);
	else to.vedges = from.vedges;
	DeepCopy(from.vcells, to.vcells, 0);
	from.enumerate_all();

	for (auto& e: to.vedges){
		for (auto& v: e->vertices) v = to.vvert[v->id];
		if (!e->left.expired()) e->left = to.vcells[e->left.lock()->id];
		if (!e->right.expired()) e->right = to.vcells[e->right.lock()->id];
	}

	for (auto& c: to.vcells)
	for (auto& e: c->edges) e = to.vedges[e->id];
}

// ===================== Miscellaneous
double HM2D::Length(const EdgeData& ed){
	double ret = 0;
	for (auto& e: ed) ret += e->length();
	return ret;
}
vector<double> HM2D::ELengths(const EdgeData& ed){
	vector<double> ret;
	for (auto& e: ed) ret.push_back(e->length());
	return ret;
}
BoundingBox HM2D::BBox(const VertexData& vd, double eps){
	return BoundingBox::Build(vd.begin(), vd.end(), eps);
}

BoundingBox HM2D::BBox(const EdgeData& ed, double eps){
	return BBox(AllVertices(ed), eps);
}
BoundingBox HM2D::BBox(const CellData& cd, double eps){
	return BBox(AllVertices(cd), eps);
}
// ================== scaling
ScaleBase HM2D::Scale01(EdgeData& ed, double a){
	auto av = AllVertices(ed);
	return HM2D::Scale01(av, a);
}
ScaleBase HM2D::Scale01(VertexData& av, double a){
	return ScaleBase::p_doscale(av.begin(), av.end(), a);
}
void HM2D::Scale(EdgeData& ed, const ScaleBase& sc){
	auto av = AllVertices(ed);
	sc.p_scale(av.begin(), av.end());
}
void HM2D::Scale(VertexData& av, const ScaleBase& sc){
	sc.p_scale(av.begin(), av.end());
}
void HM2D::Unscale(EdgeData& ed, const ScaleBase& sc){
	auto av = AllVertices(ed);
	sc.p_unscale(av.begin(), av.end());
}
void HM2D::Unscale(VertexData& av, const ScaleBase& sc){
	sc.p_unscale(av.begin(), av.end());
}

// ================== split
vector<EdgeData> HM2D::SplitData(const EdgeData& data){
	vector<vector<int>> ee = Connectivity::EdgeEdge(data);
	vector<vector<int>> sg = HMMath::Graph::SplitGraph(ee);
	vector<EdgeData> ret;
	for (auto& g: sg){
		ret.emplace_back();
		for (auto& gg: g){
			ret.back().push_back(data[gg]);
		}
	}
	return ret;
}

vector<CellData> HM2D::SplitData(const CellData& data){
	vector<vector<int>> cc = Connectivity::CellCell(data);
	vector<vector<int>> sg = HMMath::Graph::SplitGraph(cc);
	vector<CellData> ret;
	for (auto& g: sg){
		ret.emplace_back();
		for (auto& gg: g){
			ret.back().push_back(data[gg]);
		}
	}
	return ret;
}

vector<GridData> HM2D::SplitData(const GridData& data){
	vector<CellData> c = SplitData(data.vcells);
	vector<GridData> ret(c.size());
	for (int i=0; i<c.size(); ++i){
		ret[i].vcells = std::move(c[i]);
		ret[i].vedges = AllEdges(ret[i].vcells);
		ret[i].vvert = AllVertices(ret[i].vedges);
	}
	return ret;
}
