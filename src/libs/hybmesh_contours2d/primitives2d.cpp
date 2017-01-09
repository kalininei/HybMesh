#include "primitives2d.hpp"
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

// =============== contains procedures
shared_ptr<Vertex> HM2D::Contains(const VertexData& data, const Point* pnt){
	auto fnd = std::find_if(data.begin(), data.end(),
		[&pnt](const shared_ptr<Vertex>& pd){ return pd.get() == pnt; });
	if (fnd != data.end()) return *fnd;
	else return nullptr;
}
shared_ptr<Vertex> HM2D::Contains(const EdgeData& data, const Point* pnt){
	for (auto e: data){
		if (e->vertices[0].get() == pnt) return e->vertices[0];
		if (e->vertices[1].get() == pnt) return e->vertices[1];
	}
	return nullptr;
}
shared_ptr<Edge> HM2D::Contains(const EdgeData& data, const Edge* ed){
	auto fnd = std::find_if(data.begin(), data.end(),
		[&ed](const shared_ptr<Edge>& pd){ return pd.get() == ed; });
	if (fnd != data.end()) return *fnd;
	else return nullptr;
}
shared_ptr<Vertex> HM2D::Contains(const CellData& data, const Point* pnt){
	return Contains(AllVertices(data), pnt);
}
shared_ptr<Edge> HM2D::Contains(const CellData& data, const Edge* ed){
	return Contains(AllEdges(data), ed);
}
shared_ptr<Cell> HM2D::Contains(const CellData& data, const Cell* c){
	auto fnd = std::find_if(data.begin(), data.end(),
		[&c](const shared_ptr<Cell>& pd){ return pd.get() == c; });
	if (fnd != data.end()) return *fnd;
	else return nullptr;
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
std::tuple<int, double, double>
HM2D::FindClosestEdge(const EdgeData& dt, const Point& p){
	std::tuple<int, double, double> ret;
	int& ind = std::get<0>(ret);
	double& dist = std::get<1>(ret);
	double& ksi = std::get<2>(ret);
	dist = 1e99; ksi = 1e99; ind = -1;

	double k;
	for (int i=0; i<dt.size(); ++i){
		auto& e = dt[i];
		double dnew = Point::meas_section(p, *e->first(), *e->last(), k);
		//if (dnew < dist){
		if (dnew - dist < geps*geps){
			dist = dnew;
			ksi = k;
			ind = i;
			if (dist<geps*geps) break;
		}
	}

	dist = sqrt(dist);
	return ret;
}

Point HM2D::FindClosestEPoint(const EdgeData& dt, const Point& p){
	auto fec = FindClosestEdge(dt, p);
	auto e = dt[std::get<0>(fec)];
	return Point::Weigh(*e->first(), *e->last(), std::get<2>(fec));
}

std::tuple<int, double> HM2D::FindClosestNode(const VertexData& dt, const Point& p){
	if (dt.size() == 0) return std::tuple<int, double>(-1, -1);
	std::tuple<int, double> ret(0, 0);
	int& bind = std::get<0>(ret);
	double& bm = std::get<1>(ret);
	bm = Point::meas(*dt[0], p);
	for (int i=1; i<dt.size(); ++i){
		double m = Point::meas(*dt[i], p);
		if (m < bm){
			bm = m;
			bind = i;
		}
	}
	bm = sqrt(bm);
	return ret;
}


// ================== scaling
ScaleBase HM2D::Scale01(EdgeData& ed){
	auto av = AllVertices(ed);
	return HM2D::Scale01(av);
}
ScaleBase HM2D::Scale01(VertexData& av){
	return ScaleBase::p_doscale(av.begin(), av.end());
}
void HM2D::Scale(EdgeData& ed, const ScaleBase& sc){
	auto av = AllVertices(ed);
	sc.p_scale(av.begin(), av.end());
}
void HM2D::Unscale(EdgeData& ed, const ScaleBase& sc){
	auto av = AllVertices(ed);
	sc.p_unscale(av.begin(), av.end());
}
