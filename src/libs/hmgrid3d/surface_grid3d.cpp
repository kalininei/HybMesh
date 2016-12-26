#include "surface_grid3d.hpp"
#include <stack>

using namespace HMGrid3D;


ShpVector<Face> Surface::allfaces() const{
	return faces;
}
ShpVector<HMGrid3D::Edge> Surface::alledges() const{
	ShpVector<Edge> ret;
	for (auto f: faces){
		auto dt = f->alledges();
		ret.insert(ret.end(), dt.begin(), dt.end());
	}
	return aa::no_dublicates(ret);
}
ShpVector<Vertex> Surface::allvertices() const{
	ShpVector<Vertex> ret;
	for (auto f: faces){
		auto dt = f->allvertices();
		ret.insert(ret.end(), dt.begin(), dt.end());
	}
	return aa::no_dublicates(ret);
}

Surface Surface::GridSurface(HMGrid3D::GridData& g, int reversetp){
	Surface ret = GridSurface(g);
	if (reversetp == 1){
		for (auto f: ret.allfaces()) if (!f->has_right_cell()) f->reverse();
	} else if (reversetp == -1){
		for (auto f: ret.allfaces()) if (!f->has_left_cell()) f->reverse();
	}
	return ret;
}
Surface Surface::GridSurface(const HMGrid3D::GridData& g){
	Surface ret;
	for (auto f: g.vfaces){
		if (f->is_boundary()) ret.faces.push_back(f);
	}
	return ret;
}

Surface Surface::FromBoundaryType(HMGrid3D::GridData& g, int btype, int reversetp){
	Surface ret;
	for (auto f: g.vfaces){
		if (f->is_boundary() && f->boundary_type == btype)
			ret.faces.push_back(f);
	}
	if (reversetp != 0){
		if (reversetp == 1){
			for (auto f: ret.allfaces()) if (!f->has_right_cell()) f->reverse();
		} else if (reversetp == -1){
			for (auto f: ret.allfaces()) if (!f->has_left_cell()) f->reverse();
		}
	}
	return ret;
}
std::map<int, Surface> Surface::ByBoundaryTypes(HMGrid3D::GridData& g, int reversetp){
	std::map<int, Surface> ret;
	for (auto f: g.vfaces){
		if (f->is_boundary() && ret.find(f->boundary_type) == ret.end()){
			ret.emplace(f->boundary_type, FromBoundaryType(g, f->boundary_type, reversetp));
		}
	}
	return ret;
}
std::map<int, Surface> Surface::ByBoundaryTypes(const HMGrid3D::GridData& g){
	auto gp = const_cast<HMGrid3D::GridData*>(&g);
	return ByBoundaryTypes(*gp, 0);
}

Surface Surface::SubSurface(const Surface& s, Vertex* v){
	//subdivide
	vector<ShpVector<Face>> subd = Face::SubDivide(s.allfaces());
	//find part which contains v
	for (auto& fvec: subd){
		Surface ret;
		ret.faces = fvec;
		auto av = ret.allvertices();
		if (aa::shp_find(av.begin(), av.end(), v) != av.end()) return ret;
	}
	return Surface();
}

vector<Surface> Surface::AllSubSurfaces(const Surface& s){
	//subdivide
	vector<ShpVector<Face>> subd = Face::SubDivide(s.allfaces());
	//assemble surfaces
	vector<Surface> ret(subd.size());
	for (size_t i=0; i<subd.size(); ++i){
		std::swap(ret[i].faces, subd[i]);
	}
	return ret;
}

namespace{
//build table and finds edge
vector<vector<int>> face_edge_to_int(const Surface& s, const HMGrid3D::Edge* ed0, int& edcount, int& ied0){
	vector<vector<int>> ret(s.faces.size());
	auto edges = s.alledges();
	auto _indexer = aa::ptr_container_indexer(edges);
	_indexer.convert();
	int iface = 0;
	for (auto f: s.allfaces()){
		ret[iface].reserve(f->n_edges());
		for (auto e: f->edges){
			ret[iface].push_back(_indexer.index(e.get()));
		}
		++iface;
	}
	//supplementary
	edcount = edges.size();
	ied0 = aa::shp_find(edges.begin(), edges.end(), ed0) - edges.begin();
	assert(ied0 >=0 && ied0 < edcount);
	return ret;
}

//resulting face_face match face_edge order containing -1 if i-th face edge is boundary
vector<vector<int>> to_face_face(int nedges, const vector<vector<int>>& face_edge){
	vector<vector<int>> ret(face_edge.size());
	//edge face table
	vector<int> edge_face1(nedges, -1), edge_face2(nedges, -1);
	int nface = 0;
	for (auto& it1: face_edge){
		for (int e: it1){
			if (edge_face1[e]<0) edge_face1[e] = nface;
			else{
				assert(edge_face2[e]<0);
				edge_face2[e] = nface;
			}
		}
		++nface;
	}
	
	for (int iface=0; iface<face_edge.size(); ++iface){
		for (int it=0; it<face_edge[iface].size(); ++it){
			int iedge = face_edge[iface][it];
			int sib1 = edge_face1[iedge], sib2 = edge_face2[iedge];
			int sib = (sib1 == iface) ? sib2 : sib1;
			ret[iface].push_back(sib);
		}
	}

	return ret;
}

class RearrangeS{
	//temporary data
	vector<bool> used_faces;
	vector<bool> used_edges;
	vector<vector<int>> face_edge;
	vector<vector<int>> face_face;
	int _nedges;
	const Surface* input_surface;
	//resulting data
	vector<int> face_order;
	vector<int> face_edge_start;
	vector<bool> edge_rev;

	//data funcs
	int n_faces() const { return face_edge.size(); }
	int n_edges() const { return _nedges; }
	int face_size(int f) const { return face_edge[f].size(); }
	
	//algorithm functions
	void BuildProcess(int estart){
		int f1 = 0;
		for (auto& f: face_edge){
			auto fnd = std::find(f.begin(), f.end(), estart);
			if (fnd != f.end()){
				ProcessFace(f1, estart);
				break;
			}
			++f1;
		}
		assert(f1 != face_edge.size());
		//supplement face_order if not all faces were connected to estart
		if (face_order.size() != n_faces()){
			for (int i = 0; i<n_faces(); ++i) if (!used_faces[i]){
				face_order.push_back(i);
			}
			assert(face_order.size() == n_faces());
		}
	}

	bool NeedRevert(int f, int ed){
		//f - global face, ed - local edge
		auto face = input_surface->faces[f];
		auto e1 = face->edges[ed];
		auto e2 = face->edges[(ed + 1) % face->n_edges()];
		return (e1->last() == e2->first() || e1->last() == e2->last());
	}

	void ProcessFace(int f, int e){
		//f - global face, e - global edge
		if (used_faces[f]) return;
		used_faces[f] = true;
	
		//1. order
		face_order.push_back(f);
	
		//2. set start edge
		auto fnd = std::find(face_edge[f].begin(), face_edge[f].end(), e);
		assert(fnd != face_edge[f].end());
		int e0index = fnd - face_edge[f].begin();
		face_edge_start[f] = e0index;

		//3. assemble true edge order
		vector<int> loced; loced.reserve(face_edge[f].size());
		for (int i=0; i<face_size(f); ++i){
			loced.push_back((i+e0index) % face_size(f));
		}

		//4. revert possible edges
		for (int i: loced){
			int ednum = face_edge[f][i];
			if (!used_edges[ednum]){
				used_edges[ednum] = true;
				if (NeedRevert(f, i)) edge_rev[ednum] = true;
			}
		}
	
		//5. call siblings
		for (int i: loced){
			int sibface = face_face[f][i];
			if (sibface >= 0) ProcessFace(sibface, face_edge[f][i]);
		}
	}

	//apply functions
	void ReverseEdges(Surface& s){
		auto it = edge_rev.begin();
		for (auto e: s.alledges()) if (*it++) e->reverse();
	}
	void RotateFaces(Surface& s){
		auto it = face_edge_start.begin();
		for (auto f: s.allfaces()){
			std::rotate(f->edges.begin(), f->edges.begin() + *it++, f->edges.end());
		}
	}
	void PermuteFaces(Surface& s){
		ShpVector<Face> newf; newf.reserve(s.faces.size());
		for (int i: face_order) newf.push_back(s.faces[i]);
		std::swap(s.faces, newf);
	}
public:
	RearrangeS(const Surface& s, const HMGrid3D::Edge* ed0): input_surface(&s){
		int first_edge;
		//assemble connectivity tables
		face_edge = face_edge_to_int(s, ed0, _nedges, first_edge);
		face_face = to_face_face(_nedges, face_edge);
		//init tmp data
		used_faces.resize(n_faces(), false);
		used_edges.resize(n_edges(), false);
		//initial result values
		face_order.reserve(n_faces());
		face_edge_start.resize(n_faces(), 0);
		edge_rev.resize(n_edges(), false);
		//invoke recursive algorithm
		BuildProcess(first_edge);
	}

	void Apply(Surface& s){
		ReverseEdges(s);
		RotateFaces(s);
		PermuteFaces(s);
	}
};

}//anonymous namespace

void Surface::FaceRearrange(Surface& s, const Edge* ed){
	RearrangeS r(s, ed);
	r.Apply(s);
}

// ================= Algos
bool Surface::MatchTopology(const Surface& a, const Surface& b){
	auto af = a.allfaces();
	auto bf = b.allfaces();
	if (af.size() != bf.size()) return false;
	auto av = a.allvertices();
	auto bv = b.allvertices();
	if (av.size() != bv.size()) return false;

	//face->vertices connectivity
	std::vector<ShpVector<Vertex>> face_vert1, face_vert2;
	for (auto f: af) face_vert1.push_back(f->sorted_vertices());
	for (auto f: bf) face_vert2.push_back(f->sorted_vertices());

	//compare face_vert1/2
	auto _vert1 = aa::ptr_container_indexer(av); _vert1.convert();
	auto _vert2 = aa::ptr_container_indexer(bv); _vert2.convert();
	for (int i=0; i<face_vert1.size(); ++i){
		if (face_vert1[i].size() != face_vert2[i].size()) return false;
		for (int j=0; j<face_vert1[i].size(); ++j){
			int ind1 = _vert1.index(face_vert1[i][j]);
			int ind2 = _vert2.index(face_vert2[i][j]);
			if (ind1 != ind2) return false;
		}
	}
	
	return true;
}

ShpVector<HMGrid3D::Edge> Surface::ExtractBoundary(const Surface& a, Vertex v){
	auto fc = a.allfaces();
	//leave only edges which have single face connection
	ShpVector<Edge> bedges;
	for (auto& fe: Face::Connectivity::EdgeFace(fc)){
		if (fe.second.size() == 1) bedges.push_back(fe.first);
	}
	return Edge::Connect(bedges, v);
}
std::array<vector<EdgeData>, 2> Surface::ExtractAllBoundaries(const Surface& a, Vect3 right_normal){
	std::array<vector<EdgeData>, 2> ret;
	//leave only boundary edges
	EdgeData bedges;
	vector<bool> reversed;
	for (auto& fe: Face::Connectivity::EdgeFaceExtended(a.faces)){
		if (fe.second.size() > 1) continue;
		bedges.push_back(fe.first);
		reversed.push_back(!std::get<2>(fe.second[0]));
	}

	//temporary reverse edges
	auto reverse_edge = [](EdgeData& ed, vector<bool>& need){
		for (int i=0; i<ed.size(); ++i)
			if (need[i]) ed[i]->reverse();
	};
	reverse_edge(bedges, reversed);
	
	auto assemble = [](EdgeData& ed)->vector<EdgeData>{
		//enumarate vertices
		int k=0;
		for (auto& e: ed) {e->last()->id = -1; }
		for (auto& e: ed) {e->first()->id = k++; }
		vector<bool> used(ed.size(), false);
		//process
		vector<EdgeData> ret;
		while (1){
			auto fnd = std::find(used.begin(), used.end(), false);
			if (fnd == used.end()) break;
			ret.emplace_back();
			auto& rs = ret.back();

			int istart = fnd - used.begin();
			int icur = istart;
			do{
				rs.push_back(ed[icur]);
				used[icur] = true;
				icur = ed[icur]->last()->id;
			} while (istart != icur && icur>=0);

			if (icur<0) ret.resize(ret.size()-1);
		}
		return ret;
	};
	//connect all
	vector<EdgeData> alllines = assemble(bedges);

	//choose internal/external
	for (auto& line: alllines){
		assert(line.size()>0);
		auto p0 = line[0]->first();
		double area=0.;
		for (int i=1; i<line.size()-1; ++i){
			auto p1 = line[i]->first();
			auto p2 = line[i]->last();
			area += vecDot(vecCross(*p1-*p0, *p2-*p0), right_normal);
		}
		if (area > 0) ret[0].push_back(line);
		else ret[1].push_back(line);
	}
	
	//reverse edges back
	reverse_edge(bedges, reversed);

	return ret;
}

vector<Surface> Surface::ExtractSmooth(const Surface& s, double angle){
	enumerate_ids_pvec(s.faces);
	vector<Vect3> normals(s.faces.size());
	for (int i=0; i<s.faces.size(); ++i) normals[i] = s.faces[i]->left_normal();

	std::list<shared_ptr<HMGrid3D::Face>> unused(s.faces.begin(), s.faces.end());
	vector<Surface> ret;

	auto use_face = [&](std::list<shared_ptr<HMGrid3D::Face>>::iterator it, Surface* srf){
		srf->faces.push_back(*it);
		return unused.erase(it);
	};

	double badcos = cos(angle*M_PI/180.);
	while (unused.size()>0){
		ret.push_back(Surface());
		auto& news = ret.back();

		Vect3 normal1 = normals[(*unused.begin())->id];
		use_face(unused.begin(), &news);

		auto it = unused.begin();
		while (it != unused.end()){
			Vect3& normal2 = normals[(*it)->id];
			double cos = vecDot(normal1, normal2);
			if (cos > badcos){
				it = use_face(it, &news);
			} else ++it;
		}
	}

	vector<Surface> ret1;
	for (auto& r: ret){
		vector<Surface> ss = Surface::AllSubSurfaces(r);
		ret1.resize(ret1.size()+ss.size());
		for (int i=0; i<ss.size(); ++i){
			std::swap(ret1[ret1.size()-ss.size()+i].faces, ss[i].faces);
		}
	}
	
	return ret1;
}

double Surface::Volume(const FaceData& fc){
	double ret = 0.;
	Point3 p0;
	for (auto& f: fc)
	for (auto& e: f->edges) 
	for (auto& v: e->vertices){
		p0 = *v;
		goto P0FOUND;
	}
	return ret;
P0FOUND:
	for (auto& f: fc){
		auto dv = f->sorted_vertices();
		if (dv.size()<3) continue;
		for(auto p: dv) (*p) -= p0;
		Point3* p1 = dv[0].get();
		for (size_t i=1; i<dv.size()-1; ++i){
			Point3* p2 = dv[i].get();
			Point3* p3 = dv[i+1].get();
			ret += tetrahedron_volume_0(*p1, *p2, *p3);
		}
		for(auto p: dv) (*p) += p0;
	}
	return ret;
}

SurfaceTree SurfaceTree::Assemble(const Surface& idata){
	// !!! only bounding box check here
	vector<Surface> data = Surface::AllSubSurfaces(idata);
	vector<Surface*> closed;
	vector<BoundingBox3D> closed_bboxes;
	vector<Surface*> open;
	//find all closed contours, assemble bounding boxes for them
	for (int i=0; i<data.size(); ++i){
		auto ae = data[i].alledges();
		constant_ids_pvec(ae, 0);
		for (auto f: data[i].faces){
			for (auto e: f->edges){ ++e->id; }
		}
		bool isclosed = true;
		for (auto e: ae) if (e->id < 2){
			isclosed = false;
			break;
		}
		if (isclosed) closed.push_back(&data[i]);
		else open.push_back(&data[i]);

		if (isclosed){
			HMGrid3D::VertexData av;
			for (auto e: ae) constant_ids_pvec(e->vertices, 0);
			for (auto e: ae) for (auto v: e->vertices) if (v->id==0){
				av.push_back(v);
				v->id = 1;
			}
			closed_bboxes.push_back(BoundingBox3D(av));
		}
	}
	//add open contours to result
	SurfaceTree ret;
	for (auto s: open){ ret.nodes.emplace_back(new SurfaceTree::TNode(*s, -1)); }
	if (closed.size() == 0) return ret;
	//add closed contours
	for (auto s: closed) ret.nodes.emplace_back(new SurfaceTree::TNode(*s, 0));
	//build tree
	for (int i=1; i<closed.size(); ++i){
		vector<int> inners, outers;
		for (int j=0; j<i; ++j){
			int pos = closed_bboxes[i].relation(closed_bboxes[j]);
			if (pos == 1){
				inners.push_back(j);
			} else if (pos == 2){
				outers.push_back(j);
			} else if (pos == 3){
				//just a sibling, do nothing
			} else throw std::runtime_error(
				"Failed to assemble surface tree for "
				"complicated geometry");
		}
		//look up
		int lowest_outer=-1;
		int klowest = -1;
		for (int k: outers){
			int olv = ret.nodes[open.size()+k]->level;
			if (olv > lowest_outer){
				lowest_outer = olv;
				klowest = k;
			}
		}
		if (klowest >= 0){
			auto out = ret.nodes[open.size() + klowest];
			auto self = ret.nodes[open.size() + i];
			out->children.push_back(self);
			self->parent = out;
			self->level = lowest_outer+1;
		}
		//look down
		for (int k: inners){
			auto in = ret.nodes[open.size() + k];
			auto self = ret.nodes[open.size() + i];
			++(in->level);
			if (in->level == self->level + 1){
				if (!in->parent.expired()){
					auto par = in->parent.lock();
					auto fnd = std::find_if(par->children.begin(), par->children.end(),
							[&in](weak_ptr<SurfaceTree::TNode> e){
								return (e.lock() == in);
							});
					assert(fnd != par->children.end());
					par->children.erase(fnd);
				}
				in->parent = self;
			}
		}
	}
	return ret;
}

vector<SurfaceTree> SurfaceTree::crop_level1() const{
	vector<SurfaceTree> ret;
	for (auto nd: nodes) if (nd->level>=0 && nd->level % 2 == 0) {
		ret.push_back(SurfaceTree());
		auto& st = ret.back();
		st.nodes.reserve(1 + nd->children.size());
		st.nodes.emplace_back(new SurfaceTree::TNode(*nd, 0));
		auto& sr = st.nodes[0];
		for (auto ch: nd->children){
			st.nodes.emplace_back(new SurfaceTree::TNode(*ch.lock(), 1));
			auto& sch = st.nodes.back();
			sch->parent = sr;
			sr->children.push_back(sch);
		}
	}
	return ret;
}
ShpVector<SurfaceTree::TNode> SurfaceTree::roots() const {
	ShpVector<TNode> ret;
	for (int n=0; n<nodes.size(); ++n){
		if (nodes[n]->level == 0) ret.push_back(nodes[n]);
	}
	return ret;
}
HMGrid3D::EdgeData SurfaceTree::alledges() const{
	HMGrid3D::EdgeData ret;
	for (auto& n: nodes){
		auto ae = n->alledges();
		ret.resize(ret.size() + ae.size());
		std::copy(ae.begin(), ae.end(), ret.end() - ae.size());
	}
	return ret;
}
HMGrid3D::VertexData SurfaceTree::allvertices() const{
	HMGrid3D::VertexData ret;
	for (auto& n: nodes){
		auto ae = n->allvertices();
		ret.resize(ret.size() + ae.size());
		std::copy(ae.begin(), ae.end(), ret.end() - ae.size());
	}
	return ret;
}

HMGrid3D::FaceData SurfaceTree::allfaces() const{
	FaceData ret;
	for (auto& n: nodes){
		ret.insert(ret.end(), n->faces.begin(), n->faces.end());
	}
	return ret;
}

SurfTReverter::SurfTReverter(const Surface& srf): SurfTReverterBase(){
	obj = const_cast<Surface*>(&srf);
	need_revert = vector<bool>(srf.faces.size(), false);
	auto ae = srf.alledges();
	enumerate_ids_pvec(ae);
	enumerate_ids_pvec(srf.faces);
	vector<int> adjfaces1(ae.size(), -1), adjfaces2(ae.size(), -1);
	vector<int> adjlocind1(ae.size(), -1), adjlocind2(ae.size(), -1);
	for (size_t i=0; i<srf.faces.size(); ++i)
	for (size_t j=0; j<srf.faces[i]->edges.size(); ++j){
		int eid = srf.faces[i]->edges[j]->id;
		if (adjfaces1[eid] == -1){
			adjfaces1[eid] = i;
			adjlocind1[eid] = j;
		} else if (adjfaces2[eid] == -1){
			adjfaces2[eid] = i;
			adjlocind2[eid] = j;
		} else throw std::runtime_error("edge has more than two adjacent faces");
	}
	vector<bool> processed_faces(srf.faces.size(), false);
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
			auto fc = srf.faces[std::get<1>(finfo)];
			bool really_positive = fc->is_positive_edge(iedge);
			if (should_be_positive != really_positive){
				need_revert[fc->id] = true;
			}
			face_list.pop();
			for (int j=0; j<fc->n_edges(); ++j){
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
	revert();
}
void SurfTReverter::_do_revert(){
	for (size_t i=0; i<obj->faces.size(); ++i)
	if (need_revert[i]){
		obj->faces[i]->reverse();
	}
}
void SurfTReverter::_undo_revert(){
	for (size_t i=0; i<obj->faces.size(); ++i)
	if (need_revert[i]){
		obj->faces[i]->reverse();
	}
}
void SurfTReverter::reverse_all(){
	for (size_t i=0; i<need_revert.size(); ++i)
		need_revert[i] = !need_revert[i];
	if (reverted){
		for (auto& f: obj->faces) f->reverse();
	}
}
SurfTreeTReverter::SurfTreeTReverter(const SurfaceTree& srf): SurfTReverterBase(){
	obj = const_cast<SurfaceTree*>(&srf);
	for (auto& n: obj->nodes){
		if (n->isopen()) openrevs.emplace_back(new SurfTReverter(*n));
		else closedrevs.emplace_back(new SurfTReverter(*n));
	}
	int k=0;
	for (auto& n: obj->nodes) if (n->isclosed()){
		auto r = closedrevs[k++];
		double volume = Surface::Volume(*n);
		if (volume < 0 && n->level % 2 == 0) r->reverse_all();
		else if (volume > 0 && n->level % 2 == 1) r->reverse_all();
	}
	revert();
}
void SurfTreeTReverter::_do_revert(){
	for (auto r: openrevs) r->revert();
	for (auto r: closedrevs) r->revert();
}
void SurfTreeTReverter::_undo_revert(){
	for (auto r: openrevs) r->revert_back();
	for (auto r: closedrevs) r->revert_back();
}
