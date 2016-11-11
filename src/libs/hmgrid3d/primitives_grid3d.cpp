#include "primitives_grid3d.hpp"
#include "addalgo.hpp"

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
	std::map<shared_ptr<Vertex>, vector<int>> vertex_edge = Edge::Connectivity::EndVertexEdge(data);

	//4. assembling
	vector<int> assembled(1, fnd - data.begin());
	shared_ptr<Vertex> curv(v);
	while(1){
		const shared_ptr<Edge>& curedge = data[assembled.back()];
		assert(curedge->first() == curv);
		shared_ptr<Vertex> nextv = curedge->last();
		auto& ve = vertex_edge[nextv];
		assert(ve.size() == 1 || ve.size() == 2);
		int i_nextedge;
		if (ve.size() == 1) break;
		else{
			i_nextedge = ve[0];
			if (i_nextedge == assembled.back()) i_nextedge = ve[1];
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

std::map<shared_ptr<Vertex>, vector<int>>
Edge::Connectivity::EndVertexEdge(const ShpVector<Edge>& data){
	std::map<shared_ptr<Vertex>, vector<int>> ret;
	for (int i=0; i<data.size(); ++i){
		auto emp1 = ret.emplace(data[i]->first(), std::vector<int>());
		emp1.first->second.push_back(i);
		auto emp2 = ret.emplace(data[i]->last(), std::vector<int>());
		emp2.first->second.push_back(i);
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

std::map<shared_ptr<Edge>, vector<int>>
Face::Connectivity::EdgeFace(const ShpVector<Face>& data){
	std::map<shared_ptr<Edge>, vector<int>> ret;
	for (int i=0; i<data.size(); ++i){
		for (auto e: data[i]->edges){
			auto emp = ret.emplace(e, vector<int>());
			emp.first->second.push_back(i);
		}
	}
	return ret;
}

vector<vector<int>>
Face::Connectivity::FaceFace(const ShpVector<Face>& data){
	return FaceFace(EdgeFace(data), data.size());
}

vector<vector<int>>
Face::Connectivity::FaceFace(const std::map<shared_ptr<Edge>, vector<int>>& edge_face, int nfaces){
	vector<vector<int>> ret(nfaces);
	for (auto& it: edge_face){
		for (int i=0; i<it.second.size(); ++i){
			int f1 = it.second[i];
			for (int j=i+1; j<it.second.size(); ++j){
				int f2 = it.second[j];
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

void GridData::enumerate_all() const{
	enumerate_ids_pvec(vvert);
	enumerate_ids_pvec(vedges);
	enumerate_ids_pvec(vfaces);
	enumerate_ids_pvec(vcells);
}

void GridData::fill_from_cells(){
	vvert.clear();
	vfaces.clear();
	vedges.clear();

	ShpVector<Face> af;
	for (auto& c: vcells){
		std::copy(c->faces.begin(), c->faces.end(), std::back_inserter(af));
	}
	ShpVector<Edge> ae;
	for (auto& f: af){
		std::copy(f->edges.begin(), f->edges.end(), std::back_inserter(ae));
	}
	ShpVector<Vertex> av;
	for (auto& e: ae){
		std::copy(e->vertices.begin(), e->vertices.end(), std::back_inserter(av));
	}

	constant_ids_pvec(af, -1);
	constant_ids_pvec(ae, -1);
	constant_ids_pvec(av, -1);

	for (auto f: af) if (f->id == -1){
		vfaces.push_back(f);
		f->id = 0;
	}
	for (auto e: ae) if (e->id == -1){
		vedges.push_back(e);
		e->id = 0;
	}
	for (auto v: av) if (v->id == -1){
		vvert.push_back(v);
		v->id = 0;
	}
}

void GridData::add_alldata(const GridData& from){
	_DUMMY_FUN_;
}

//reallocation
void HMGrid3D::ReallocateAll(VertexData& vd){
	for (auto& v: vd) v.reset(new Vertex(*v));
}

void HMGrid3D::ReallocateAll(EdgeData& ed){
	for (auto e: ed)
	for (auto v: e->vertices) v->id = -2;
	VertexData uvert;
	for (auto e: ed)
	for (auto v: e->vertices) if (v->id == -2){
		uvert.push_back(v);
		v->id = -1;
	}
	VertexData cpvert = uvert;
	ReallocateAll(cpvert);

	enumerate_ids_pvec(uvert);
	for (auto& e: ed){
		shared_ptr<HMGrid3D::Edge> ne(new Edge(*e));
		for (int i=0; i<e->vertices.size(); ++i){
			ne->vertices[i] = cpvert[ne->vertices[i]->id];
		}
		e = ne;
	}
}
void HMGrid3D::ReallocateAll(FaceData& fd){
	//unique edges and nodes
	for (auto f: fd)
	for (auto e: f->edges) e->id = -2;

	EdgeData uedges;

	for (auto f: fd)
	for (auto e: f->edges) if (e->id == -2){
		uedges.push_back(e);
		e->id = -1;
	}
	
	EdgeData cpedges = uedges;
	ReallocateAll(cpedges);

	enumerate_ids_pvec(uedges);
	for (auto& f: fd){
		shared_ptr<HMGrid3D::Face> nf(new Face(*f));
		for (int i=0; i<nf->n_edges(); ++i){
			nf->edges[i] = cpedges[nf->edges[i]->id];
		}
		f = nf;
	}
}
