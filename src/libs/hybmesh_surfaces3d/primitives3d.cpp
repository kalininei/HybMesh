#include "primitives3d.hpp"
#include "addalgo.hpp"
#include "debug3d.hpp"
#include "contabs3d.hpp"
#include "surface.hpp"

using namespace HM3D;

std::tuple<Vertex*, int, double>
HM3D::FindClosestVertex(const VertexData& vec, Vertex v){
	std::tuple<Vertex*, int, double> ret;
	Vertex*& vout = std::get<0>(ret);
	int& ind = std::get<1>(ret);
	double& meas = std::get<2>(ret);

	vout = 0; ind = 0; meas = -1;
	int it=0;
	for (auto& x: vec){
		double m = Vertex::meas(v, *x);
		if (vout == 0 || m<meas){
			vout = x.get(); meas = m; ind = it;
		}
		++it;
	}

	return ret;
}

// ================= Edge
EdgeData HM3D::Connect(const EdgeData& data, Vertex app_v){
	//1. find closest vertex to v
	ShpVector<Vertex> vall;
	for (auto d: data) { vall.push_back(d->first()); vall.push_back(d->last()); }
	vall = aa::no_dublicates(vall);
	auto fres = FindClosestVertex(vall, app_v);
	const shared_ptr<Vertex> v = vall[std::get<1>(fres)];

	//2. find edge which includes v as start node
	auto fnd = std::find_if(data.begin(), data.end(),
			[&v](shared_ptr<Edge> d){ return d->first() == v;});
	assert(fnd != data.end());

	//3. assemble vertex_edge connectivity
	auto vertex_edge = Connectivity::VertexEdge(data);
	for (int i=0; i<vertex_edge.size(); ++i) vertex_edge[i].v->id = i;

	//4. assembling
	vector<int> assembled(1, fnd - data.begin());
	shared_ptr<Vertex> curv(v);
	while(1){
		const shared_ptr<Edge>& curedge = data[assembled.back()];
		assert(curedge->first() == curv);
		shared_ptr<Vertex> nextv = curedge->last();
		auto& ve = vertex_edge[nextv->id];
		assert(ve.size() == 1 || ve.size() == 2);
		int i_nextedge;
		if (ve.size() == 1) break;
		else{
			i_nextedge = ve.eind[0];
			if (i_nextedge == assembled.back()) i_nextedge = ve.eind[1];
		}
		if (i_nextedge == assembled[0]) break;
		else assembled.push_back(i_nextedge);
		std::swap(curv, nextv);
	}

	//5. write return vector
	ShpVector<Edge> ret; ret.reserve(assembled.size());
	for (int i: assembled) ret.push_back(data[i]);
	return ret;
}

double Edge::measure() const{
	return Vertex::meas(*first(), *last());
}
double Edge::length() const{
	return sqrt(measure());
}
void Edge::reverse(){
	std::reverse(vertices.begin(), vertices.end());
}

// ================= Face
ShpVector<Vertex> Face::sorted_vertices() const{
	assert(edges.size() > 1);
	ShpVector<Vertex> ret;
	auto it = edges.begin();
	auto it2 = std::next(it);
	while (it2 != edges.end()){
		if ((*it)->vertices.back() == (*it2)->vertices[0] ||
		    (*it)->vertices.back() == (*it2)->vertices.back()){
			ret.insert(ret.end(), (*it)->vertices.begin(), (*it)->vertices.end()-1);
		} else {
			ret.insert(ret.end(), (*it)->vertices.rbegin(), (*it)->vertices.rend()-1);
		}
		++it; ++it2;
	}
	it2 = edges.begin();
	if ((*it)->vertices.back() == (*it2)->vertices[0] ||
	    (*it)->vertices.back() == (*it2)->vertices.back()){
		ret.insert(ret.end(), (*it)->vertices.begin(), (*it)->vertices.end()-1);
	} else {
		ret.insert(ret.end(), (*it)->vertices.rbegin(), (*it)->vertices.rend()-1);
	}
	return ret;
}

bool Face::change_edge(shared_ptr<Edge> from, shared_ptr<Edge> to){
	int ind = std::find(edges.begin(), edges.end(), from) - edges.begin();
	if (ind == edges.size()) return false;
	//underlying vertices
	auto from1 = from->first();
	auto from2 = from->last();
	auto to1 = to->first();
	auto to2 = to->last();
	if (*from1 != *to1) std::swap(to1, to2);
	assert(*from1 == *to1 && *from2 == *to2);

	edges[ind] = to;
	auto edgeprev = edges[ (ind==0)?edges.size()-1:ind-1 ];
	auto edgenext = edges[ (ind==edges.size()-1)?0:ind+1 ];

	if (edgeprev->first() == from1) edgeprev->vertices[0] = to1;
	else if (edgeprev->last() == from1) edgeprev->vertices.back() = to1;
	else if (edgeprev->first() == from2) edgeprev->vertices[0] = to2;
	else if (edgeprev->last() == from2) edgeprev->vertices.back() = to2;

	if (edgenext->first() == from1) edgenext->vertices[0] = to1;
	else if (edgenext->last() == from1) edgenext->vertices.back() = to1;
	else if (edgenext->first() == from2) edgenext->vertices[0] = to2;
	else if (edgenext->last() == from2) edgenext->vertices.back() = to2;

	return true;
}


void Face::reverse(){
	std::swap(left, right);
	std::reverse(edges.begin(), edges.end());
}

void Face::correct_edge_directions(){
	assert(edges.size()>1);
	//first edge
	auto e1 = edges.back();
	auto e2 = edges[0];
	if (e2->first() != e1->last() && e2->first() != e1->first()) e2->reverse();
	//other edges
	for (auto i=1; i<edges.size(); ++i){
		auto e1 = edges[i-1], e2 = edges[i];
		if (e2->first() != e1->last()) e2->reverse();
	}
}
bool Face::is_positive_edge(int eindex){
	int enext = eindex + 1;
	if (enext == edges.size()) enext = 0;
	if (edges[eindex]->last() == edges[enext]->first()) return true;
	if (edges[eindex]->last() == edges[enext]->last()) return true;
	return false;
}
std::array<Point3, 3> Face::mean_points() const{
	//quick procedure
	assert(edges.size() > 2);
	auto op = sorted_vertices();
	if (edges.size() == 3){
		return {*op[0], *op[1], *op[2]};
	} else if (edges.size() == 4){
		return {*op[0], *op[1], (*op[2] + *op[3])/2};
	} else {
		for (int i=2; i<op.size(); ++i){
			double cp = vecLen(vecCross(*op[1]-*op[0], *op[i]-*op[0]));
			if (!ISZERO(cp)){
				double sin = cp/vecLen(*op[1]-*op[0])/vecLen(*op[i]-*op[0]);
				if (sin > 0){
					return {*op[0], *op[1], *op[i]};
				} else if (sin<0){
					return {*op[0], *op[i], *op[1]};
				}
			}
		}
	}
	throw std::runtime_error("Can not compute mean plane of the face");
}
Vect3 Face::left_normal() const{
	std::array<Point3, 3> mp = mean_points();
	auto ret = vecCross(mp[2]-mp[0], mp[1]-mp[0]);
	vecNormalize(ret);
	return ret;
}

// ================== Cell
std::tuple<int, int, int> Cell::n_fev() const{
	int nf = faces.size();
	int ne = 0;
	int nv = 0;
	for (auto& f: faces)
	for (auto& e: f->edges){
		e->id=0;
		for (auto& v: e->vertices) v->id=0;
	}
	for (auto& f: faces)
	for (auto& e: f->edges) if (e->id == 0){
		e->id=1;
		++ne;
		for (auto& v: e->vertices) if (v->id==0){
			v->id=1;
			++nv;
		}
	}
	return std::make_tuple(nf, ne, nv);
}
bool Cell::change_face(shared_ptr<Face> from, shared_ptr<Face> to){
	int ind = std::find(faces.begin(), faces.end(), from) - faces.begin();
	if (ind == faces.size()) return false;

	//underlying edges
	std::vector<shared_ptr<Edge>> from_edges(from->edges.begin(), from->edges.end());
	std::vector<shared_ptr<Edge>> to_edges(to->edges.begin(), to->edges.end());
	auto& p1 = *(*from_edges.begin())->first();
	auto& p2 = *(*from_edges.begin())->last();
	for (auto it = to_edges.begin(); it!=to_edges.end(); ++it){
		auto& p3 = *(*it)->first();
		auto& p4 = *(*it)->last();
		if ( (p1 == p3 && p2 == p4) || (p1 == p4 && p2 == p3) ){
			std::rotate(to_edges.begin(), it, to_edges.end());
			break;
		}
	}
	//underlying vertices
	Face tmpface;
	tmpface.edges = from_edges;
	std::vector<shared_ptr<Vertex>> from_vert = tmpface.sorted_vertices();
	tmpface.edges = to_edges;
	std::vector<shared_ptr<Vertex>> to_vert = tmpface.sorted_vertices();

	for (auto f: faces)
	for (auto e: f->edges){
		e->id = -1;
		for (auto v: e->vertices) v->id=-1;
	}
	for (int i=0; i<from_edges.size(); ++i) from_edges[i]->id = i;
	for (int i=0; i<from_vert.size(); ++i) from_vert[i]->id = i;
	for (int i=0; i<to_edges.size(); ++i) to_edges[i]->id = -2;
	for (int i=0; i<to_vert.size(); ++i) to_vert[i]->id = -2;

	//change target face
	faces[ind] = to;

	for (int i=0; i<faces.size(); ++i) if (i != ind)
	for (int j=0; j<faces[i]->edges.size(); ++j){
		//change edges of non-target faces
		int& id = faces[i]->edges[j]->id;
		if (id >= 0){
			faces[i]->edges[j] =to_edges[id];
		}
		//change vertices of non-target edges
		if (id == -1){
			for (auto& v: faces[i]->edges[j]->vertices){
				int& id2 = v->id;
				if (id2 >= 0) v = to_vert[id2];
			}
			id = -2;
		}
	}

	return true;
}

double Cell::volume() const{
	double ret;

	std::vector<bool> reverted(faces.size(), false);
	for (int i=0; i<faces.size(); ++i)
	if (faces[i]->left.lock().get() != this) reverted[i] = true;
	auto facerev = [&](){
		for (int i=0; i<faces.size(); ++i)
		if (reverted[i]) faces[i]->reverse();
	};

	facerev();
	try{
		ret = HM3D::Surface::Volume(faces);
	} catch (...){
		facerev();
		throw;
	}
	facerev();

	return ret;
}

vector<double> HM3D::Volumes(const CellData& cd){
	vector<double> ret(cd.size());
	for (int i=0; i<cd.size(); ++i) ret[i] = cd[i]->volume();
	return ret;
}
double HM3D::SumVolumes(const CellData& cd){
	double ret = 0;
	for (auto x: Volumes(cd)) ret += x;
	return ret;
}

void GridData::enumerate_all() const{
	aa::enumerate_ids_pvec(vvert);
	aa::enumerate_ids_pvec(vedges);
	aa::enumerate_ids_pvec(vfaces);
	aa::enumerate_ids_pvec(vcells);
}

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

//deep copy procedures
void HM3D::DeepCopy(const VertexData& from, VertexData& to){
	shared_vec_deepcopy(from, to);
}
void HM3D::DeepCopy(const EdgeData& from, EdgeData& to, int level){
	shared_vec_deepcopy(from, to);
	if (level == 1){
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
void HM3D::DeepCopy(const FaceData& from, FaceData& to, int level){
	shared_vec_deepcopy(from, to);
	if (level > 0){
		EdgeData eorig = AllEdges(from);
		EdgeData enew;
		DeepCopy(eorig, enew, level - 1);
		aa::enumerate_ids_pvec(eorig);
		aa::enumerate_ids_pvec(enew);
		auto fit = to.end() - from.size();
		while (fit != to.end()){
			for (auto& e: (*fit++)->edges){
				e = enew[e->id];
			}
		}
	}
}
void HM3D::DeepCopy(const CellData& from, CellData& to, int level){
	shared_vec_deepcopy(from, to);
	if (level > 0){
		FaceData forig = AllFaces(from);
		FaceData fnew;
		DeepCopy(forig, fnew, level - 1);
		aa::enumerate_ids_pvec(forig);
		aa::enumerate_ids_pvec(fnew);
		auto cit = to.end() - from.size();
		int k = 0;
		while (cit != to.end()){
			for (auto& f: (*cit)->faces){
				f = fnew[f->id];
				if (forig[f->id]->left.lock() == from[k]) f->left = *cit;
				else if (forig[f->id]->right.lock() == from[k]) f->right = *cit;
			}
			++k;
			++cit;
		}
	}
}
VertexData HM3D::AllVertices(const EdgeData& from){
	for (auto& e: from)
	for (auto& v: e->vertices) v->id = 0;
	VertexData ret;
	for (auto& e: from)
	for (auto& v: e->vertices) if (v->id == 0){
		ret.push_back(v);
		v->id = 1;
	}
	return ret;
}
EdgeData HM3D::AllEdges(const FaceData& from){
	for (auto& f: from)
	for (auto& e: f->edges) e->id = 0;
	EdgeData ret;
	for (auto& f: from)
	for (auto& e: f->edges) if (e->id == 0){
		ret.push_back(e);
		e->id = 1;
	}
	return ret;
}
FaceData HM3D::AllFaces(const CellData& from){
	for (auto& c: from)
	for (auto& f: c->faces) f->id = 0;
	FaceData ret;
	for (auto& c: from)
	for (auto& f: c->faces) if (f->id == 0){
		ret.push_back(f);
		f->id = 1;
	}
	return ret;
}
VertexData HM3D::AllVertices(const FaceData& from){
	return AllVertices(AllEdges(from));
}
VertexData HM3D::AllVertices(const CellData& from){
	return AllVertices(AllEdges(AllFaces(from)));
}
EdgeData HM3D::AllEdges(const CellData& from){
	return AllEdges(AllFaces(from));
}
std::tuple<VertexData> HM3D::AllPrimitives(const EdgeData& from){
	std::tuple<VertexData> ret;
	std::get<0>(ret) = AllVertices(from);
	return ret;
}
std::tuple<VertexData, EdgeData> HM3D::AllPrimitives(const FaceData& from){
	std::tuple<VertexData, EdgeData> ret;
	for (auto& f: from)
	for (auto& e: f->edges){
		e->id = 0;
		for (auto& v: e->vertices) v->id = 0;
	}
	for (auto& f: from)
	for (auto& e: f->edges) if (e->id == 0){
		std::get<1>(ret).push_back(e);
		e->id = 1;
		for (auto& v: e->vertices) if (v->id == 0){
			std::get<0>(ret).push_back(v);
			v->id = 1;
		}
	}
	return ret;
}
std::tuple<VertexData, EdgeData, FaceData> HM3D::AllPrimitives(const CellData& from){
	std::tuple<VertexData, EdgeData, FaceData> ret;
	for (auto& c: from)
	for (auto& f: c->faces){
		f->id = 0;
		for (auto& e: f->edges){
			e->id = 0;
			for (auto& v: e->vertices) v->id = 0;
		}
	}
	for (auto& c: from)
	for (auto& f: c->faces) if (f->id == 0){
		f->id = 1;
		std::get<2>(ret).push_back(f);
		for (auto& e: f->edges) if (e->id == 0){
			std::get<1>(ret).push_back(e);
			e->id = 1;
			for (auto& v: e->vertices) if (v->id == 0){
				std::get<0>(ret).push_back(v);
				v->id = 1;
			}
		}
	}
	return ret;
}

void HM3D::DeepCopy(const GridData& from, GridData& to, int level){
	to.clear();
	if (level>=3) DeepCopy(from.vvert, to.vvert);
	else to.vvert = from.vvert;
	if (level>=2) DeepCopy(from.vedges, to.vedges, 0);
	else to.vedges = from.vedges;
	if (level>=1) DeepCopy(from.vfaces, to.vfaces, 0);
	else to.vfaces= from.vfaces;
	DeepCopy(from.vcells, to.vcells, 0);
	from.enumerate_all();

	for (auto& e: to.vedges)
	for (auto& v: e->vertices) v = to.vvert[v->id];

	for (auto& f: to.vfaces){
		for (auto& e: f->edges) e = to.vedges[e->id];
		if (!f->left.expired()) f->left = to.vcells[f->left.lock()->id];
		if (!f->right.expired()) f->right = to.vcells[f->right.lock()->id];
	}

	for (auto& c: to.vcells)
	for (auto& f: c->faces) f = to.vfaces[f->id];
}

void HM3D::Unscale2D(VertexData& vd, const ScaleBase& sc){
	for (auto& v: vd){
		v->x = (v->x * sc.L) + sc.p0.x;
		v->y = (v->y * sc.L) + sc.p0.y;
		v->z = (v->z * sc.L);
	}
}
