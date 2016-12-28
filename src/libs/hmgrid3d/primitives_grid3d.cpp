#include "primitives_grid3d.hpp"
#include "addalgo.hpp"
#include "surface_grid3d.hpp"
#include "debug_grid3d.hpp"

using namespace HMGrid3D;

double Vertex::measure(const Vertex& a, const Vertex& b){
	double xd = a.x - b.x;
	double yd = a.y - b.y;
	double zd = a.z - b.z;
	return xd*xd + yd*yd + zd*zd;
}

std::tuple<Vertex*, int, double>
Vertex::FindClosestVertex(const ShpVector<Vertex>& vec, Vertex v){
	std::tuple<Vertex*, int, double> ret;
	Vertex*& vout = std::get<0>(ret);
	int& ind = std::get<1>(ret);
	double& meas = std::get<2>(ret);

	vout = 0; ind = 0; meas = -1;
	int it=0;
	for (auto& x: vec){
		double m = Vertex::measure(v, *x);
		if (vout == 0 || m<meas){
			vout = x.get(); meas = m; ind = it;
		}
		++it;
	}

	return ret;
}

void Vertex::Unscale2D(ShpVector<Vertex>& vec, ScaleBase sc){
	for(auto v: vec){
		v->x *= sc.L; v->x += sc.p0.x;
		v->y *= sc.L; v->y += sc.p0.y;
		v->z *= sc.L;

	}
}

// ================= Edge
ShpVector<Edge> Edge::Connect(const ShpVector<Edge>& data, Vertex app_v){
	//1. find closest vertex to v
	ShpVector<Vertex> vall;
	for (auto d: data) { vall.push_back(d->first()); vall.push_back(d->last()); }
	vall = aa::no_dublicates(vall);
	auto fres = Vertex::FindClosestVertex(vall, app_v);
	const shared_ptr<Vertex> v = vall[std::get<1>(fres)];

	//2. find edge which includes v as start node
	auto fnd = std::find_if(data.begin(), data.end(),
			[&v](shared_ptr<Edge> d){ return d->first() == v;});
	assert(fnd != data.end());

	//3. assemble vertex_edge connectivity
	auto vertex_edge = Edge::Connectivity::EndVertexEdge(data);
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

vector<Edge::Connectivity::EndVertexEdgeR>
Edge::Connectivity::EndVertexEdge(const EdgeData& data){
	auto ev = AllEndVertices(data);
	vector<Edge::Connectivity::EndVertexEdgeR> ret(ev.size());
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

double Edge::measure() const{
	double ret = 0;
	for (int i = 0; i<(int)vertices.size()-1; ++i){
		ret += Vertex::measure(*vertices[i], *vertices[i+1]);
	}
	return ret;
}
double Edge::length() const{
	//!!!! this is wrong for multiple vertex edge
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


ShpVector<Vertex> Face::allvertices() const{
	ShpVector<Vertex> ret;
	for (auto e: edges){
		ret.insert(ret.end(), e->vertices.begin(), e->vertices.end());
	}
	return aa::no_dublicates(ret);
}

ShpVector<Edge> Face::alledges() const{ return edges; }

void Face::reverse(){
	std::swap(left, right);
	std::reverse(edges.begin(), edges.end());
}

void Face::correct_edge_directions(){
	assert(n_edges()>1);
	//first edge
	auto e1 = edges.back();
	auto e2 = edges[0];
	if (e2->first() != e1->last() && e2->first() != e1->first()) e2->reverse();
	//other edges
	for (auto i=1; i<n_edges(); ++i){
		auto e1 = edges[i-1], e2 = edges[i];
		if (e2->first() != e1->last()) e2->reverse();
	}
}
bool Face::is_positive_edge(int eindex){
	int enext = eindex + 1;
	if (enext == n_edges()) enext = 0;
	if (edges[eindex]->last() == edges[enext]->first()) return true;
	if (edges[eindex]->last() == edges[enext]->last()) return true;
	return false;
}
std::array<Point3, 3> Face::mean_points() const{
	//quick procedure
	assert(n_edges() > 2);
	auto op = sorted_vertices();
	if (n_edges() == 3){
		return {*op[0], *op[1], *op[2]};
	} else if (n_edges() == 4){
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

vector<Face::Connectivity::EdgeFaceR>
Face::Connectivity::EdgeFace(const FaceData& data){
	auto ae = AllEdges(data);
	vector<Face::Connectivity::EdgeFaceR> ret(ae.size());
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

vector<Face::Connectivity::EdgeFaceExtendedR>
Face::Connectivity::EdgeFaceExtended(const FaceData& data){
	auto ae = AllEdges(data);
	vector<Face::Connectivity::EdgeFaceExtendedR> ret(ae.size());
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

vector<vector<int>>
Face::Connectivity::FaceFace(const ShpVector<Face>& data){
	return FaceFace(EdgeFace(data), data.size());
}

vector<vector<int>>
Face::Connectivity::FaceFace(const vector<Face::Connectivity::EdgeFaceR>& edge_face, int nfaces){
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

namespace{
void connect_faces(int fc, const vector<vector<int>>& ff, vector<bool>& used, vector<int>& ret){
	if (used[fc]) return;
	used[fc] = true;
	ret.push_back(fc);
	for (int sibf: ff[fc]) connect_faces(sibf, ff, used, ret);
}
}
std::vector<ShpVector<Face>> Face::SubDivide(const ShpVector<Face>& fvec){
	vector<vector<int>> face_face = Face::Connectivity::FaceFace(fvec);
	vector<bool> used_faces(fvec.size(), false);
	vector<ShpVector<Face>> ret;
	while (1){
		auto fnd = std::find(used_faces.begin(), used_faces.end(), false);
		if (fnd == used_faces.end()) break;
		vector<int> fc;
		connect_faces(fnd-used_faces.begin(), face_face, used_faces, fc);
		ret.push_back(ShpVector<Face>());
		ret.back().reserve(fc.size());
		for (int i: fc) ret.back().push_back(fvec[i]);
	}
	return ret;
}
void Face::SetBoundaryTypes(const ShpVector<Face>& fvec, std::function<int(Vertex, int)> bfun){
	for (auto f: fvec) if (f->is_boundary()){
		Vertex cv;
		auto av = f->allvertices();
		for (auto v: av){ cv.x += v->x; cv.y += v->y; cv.z += v->z; }
		cv.x /= av.size(); cv.y /= av.size(); cv.z /= av.size();
		f->boundary_type = bfun(cv, f->boundary_type);
	}
}

// ================== Cell
int Cell::n_faces() const{
	return faces.size();
}

int Cell::n_edges() const{
	return alledges().size();
}
int Cell::n_vertices() const{
	return allvertices().size();
}

std::tuple<int, int, int> Cell::n_fev() const{
	std::unordered_set<shared_ptr<Vertex>> ord;
	int e2 = 0;
	for (auto f: faces){
		e2 += f->n_edges();
		for (auto e: f->edges) ord.insert(e->vertices.begin(), e->vertices.end());
	}

	return std::make_tuple(n_faces(), (int)e2/2, (int)ord.size());
}

ShpVector<Vertex> Cell::allvertices() const{
	ShpVector<Vertex> ret;
	for (auto f: faces){
		auto dt = f->allvertices();
		ret.insert(ret.end(), dt.begin(), dt.end());
	}
	return aa::no_dublicates(ret);
}

ShpVector<Face> Cell::allfaces() const{ return faces; }

ShpVector<Edge> Cell::alledges() const{
	ShpVector<Edge> ret;
	for (auto c: faces){
		auto dt = c->edges;
		ret.insert(ret.end(), dt.begin(), dt.end());
	}
	return aa::no_dublicates(ret);
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

	std::vector<bool> reverted(n_faces(), false);
	for (int i=0; i<n_faces(); ++i)
	if (faces[i]->left.lock().get() != this) reverted[i] = true;
	auto facerev = [&](){
		for (int i=0; i<faces.size(); ++i)
		if (reverted[i]) faces[i]->reverse();
	};

	facerev();
	try{
		ret = HMGrid3D::Surface::Volume(faces);
	} catch (...){
		facerev();
		throw;
	}
	facerev();

	return ret;
}

vector<double> Cell::Volumes(const CellData& cd){
	vector<double> ret(cd.size());
	for (int i=0; i<cd.size(); ++i) ret[i] = cd[i]->volume();
	return ret;
}
double Cell::SumVolumes(const CellData& cd){
	double ret = 0;
	for (auto x: Volumes(cd)) ret += x;
	return ret;
}

void GridData::enumerate_all() const{
	enumerate_ids_pvec(vvert);
	enumerate_ids_pvec(vedges);
	enumerate_ids_pvec(vfaces);
	enumerate_ids_pvec(vcells);
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
void HMGrid3D::DeepCopy(const VertexData& from, VertexData& to){
	shared_vec_deepcopy(from, to);
}
void HMGrid3D::DeepCopy(const EdgeData& from, EdgeData& to, int level){
	shared_vec_deepcopy(from, to);
	if (level == 1){
		VertexData vorig = AllVertices(from);
		VertexData vnew;
		DeepCopy(vorig, vnew);
		enumerate_ids_pvec(vorig);
		enumerate_ids_pvec(vnew);
		auto eit = to.end() - from.size();
		while (eit != to.end()){
			for (auto& v: (*eit++)->vertices){
				v = vnew[v->id];
			}
		}
	}
}
void HMGrid3D::DeepCopy(const FaceData& from, FaceData& to, int level){
	shared_vec_deepcopy(from, to);
	if (level > 0){
		EdgeData eorig = AllEdges(from);
		EdgeData enew;
		DeepCopy(eorig, enew, level - 1);
		enumerate_ids_pvec(eorig);
		enumerate_ids_pvec(enew);
		auto fit = to.end() - from.size();
		while (fit != to.end()){
			for (auto& e: (*fit++)->edges){
				e = enew[e->id];
			}
		}
	}
}
void HMGrid3D::DeepCopy(const CellData& from, CellData& to, int level){
	shared_vec_deepcopy(from, to);
	if (level > 0){
		FaceData forig = AllFaces(from);
		FaceData fnew;
		DeepCopy(forig, fnew, level - 1);
		enumerate_ids_pvec(forig);
		enumerate_ids_pvec(fnew);
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
//extract procedures
VertexData HMGrid3D::AllVertices(const EdgeData& from){
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
VertexData HMGrid3D::AllEndVertices(const EdgeData& from){
	for (auto& e: from){
		e->first()->id = 0;
		e->last()->id = 0;
	}
	VertexData ret;
	for (auto& e: from){
		if (e->first()->id == 0){
			ret.push_back(e->first());
			e->first()->id = 1;
		}
		if (e->last()->id == 0){
			ret.push_back(e->last());
			e->last()->id = 1;
		}
	}
	return ret;
}
EdgeData HMGrid3D::AllEdges(const FaceData& from){
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
FaceData HMGrid3D::AllFaces(const CellData& from){
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
VertexData HMGrid3D::AllVertices(const FaceData& from){
	return AllVertices(AllEdges(from));
}
VertexData HMGrid3D::AllEndVertices(const FaceData& from){
	return AllEndVertices(AllEdges(from));
}
VertexData HMGrid3D::AllVertices(const CellData& from){
	return AllVertices(AllEdges(AllFaces(from)));
}
VertexData HMGrid3D::AllEndVertices(const CellData& from){
	return AllEndVertices(AllEdges(AllFaces(from)));
}
EdgeData HMGrid3D::AllEdges(const CellData& from){
	return AllEdges(AllFaces(from));
}
std::tuple<VertexData> HMGrid3D::AllPrimitives(const EdgeData& from){
	std::tuple<VertexData> ret;
	std::get<0>(ret) = AllVertices(from);
	return ret;
}
std::tuple<VertexData, EdgeData> HMGrid3D::AllPrimitives(const FaceData& from){
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
std::tuple<VertexData, EdgeData, FaceData> HMGrid3D::AllPrimitives(const CellData& from){
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
