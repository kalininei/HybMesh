#include "surface_grid3d.hpp"

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

Surface Surface::FromBoundaryType(HMGrid3D::SGrid& g, int btype, int reversetp){
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
std::map<int, Surface> Surface::ByBoundaryTypes(HMGrid3D::SGrid& g, int reversetp){
	std::map<int, Surface> ret;
	for (auto f: g.vfaces){
		if (f->is_boundary() && ret.find(f->boundary_type) == ret.end()){
			ret.emplace(f->boundary_type, FromBoundaryType(g, f->boundary_type, reversetp));
		}
	}
	return ret;
}
std::map<int, Surface> Surface::ByBoundaryTypes(const HMGrid3D::SGrid& g){
	auto gp = const_cast<HMGrid3D::SGrid*>(&g);
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

namespace{
//build table and finds edge
vector<vector<int>> face_edge_to_int(const Surface& s, const HMGrid3D::Edge* ed0, int& edcount, int& ied0){
	vector<vector<int>> ret(s.n_faces());
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
		ShpVector<Face> newf; newf.reserve(s.n_faces());
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
