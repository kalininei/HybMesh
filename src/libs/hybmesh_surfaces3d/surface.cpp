#include <stack>
#include "surface.hpp"
#include "debug3d.hpp"
#include "contabs3d.hpp"

namespace hs = HM3D::Surface;
using namespace hs;
using namespace HM3D;

bool hs::IsClosed(const FaceData& fd){
	auto ae = AllEdges(fd);
	aa::constant_ids_pvec(ae, 0);
	for (auto f: fd){
		for (auto e: f->edges){ ++e->id; }
	}
	bool isclosed = true;
	for (auto e: ae) if (e->id != 2){
		isclosed = false;
		break;
	}
	return isclosed;
}

FaceData hs::GridSurface(const GridData& g){
	FaceData ret;
	for (auto f: g.vfaces){
		if (f->is_boundary()) ret.push_back(f);
	}
	return ret;
}

std::map<int, FaceData> hs::GridSurfaceBType(const HM3D::GridData& g){
	FaceData all = GridSurface(g);
	std::map<int, FaceData> ret;
	for (auto f: all){
		auto er = ret.emplace(f->boundary_type, FaceData());
		er.first->second.push_back(f);
	}
	return ret;
}

void hs::SetBoundaryTypes(const FaceData& fvec, std::function<int(Vertex, int)> bfun){
	for (auto f: fvec) if (f->is_boundary()){
		Vertex cv;
		auto av = AllVertices(f->edges);
		for (auto v: av){ cv.x += v->x; cv.y += v->y; cv.z += v->z; }
		cv.x /= av.size(); cv.y /= av.size(); cv.z /= av.size();
		f->boundary_type = bfun(cv, f->boundary_type);
	}
}


namespace{
void connect_faces(int fc, const vector<vector<int>>& ff, vector<bool>& used, vector<int>& ret){
	if (used[fc]) return;
	used[fc] = true;
	ret.push_back(fc);
	for (int sibf: ff[fc]) connect_faces(sibf, ff, used, ret);
}
}

vector<FaceData> hs::AllSubSurfaces(const FaceData& s){
	vector<vector<int>> face_face = Connectivity::FaceFace(s);
	vector<bool> used_faces(s.size(), false);
	vector<ShpVector<Face>> ret;
	while (1){
		auto fnd = std::find(used_faces.begin(), used_faces.end(), false);
		if (fnd == used_faces.end()) break;
		vector<int> fc;
		connect_faces(fnd-used_faces.begin(), face_face, used_faces, fc);
		ret.push_back(ShpVector<Face>());
		ret.back().reserve(fc.size());
		for (int i: fc) ret.back().push_back(s[i]);
	}
	return ret;
}

FaceData hs::SubSurface(const FaceData& s, const Vertex* v){
	//subdivide
	vector<FaceData> subd = AllSubSurfaces(s);
	//find part which contains v
	for (auto& fvec: subd){
		auto av = AllVertices(fvec);
		if (aa::shp_find(av.begin(), av.end(), v) != av.end()) return fvec;
	}
	return FaceData();
}


namespace{
//build table and finds edge
vector<vector<int>> face_edge_to_int(const FaceData& s, const Edge* ed0, int& edcount, int& ied0){
	vector<vector<int>> ret(s.size());
	auto edges = AllEdges(s);
	auto _indexer = aa::ptr_container_indexer(edges);
	_indexer.convert();
	int iface = 0;
	for (auto f: s){
		ret[iface].reserve(f->edges.size());
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
	const FaceData* input_surface;
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
		auto face = (*input_surface)[f];
		auto e1 = face->edges[ed];
		auto e2 = face->edges[(ed + 1) % face->edges.size()];
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
	void ReverseEdges(FaceData& s){
		auto it = edge_rev.begin();
		for (auto e: AllEdges(s)) if (*it++) e->reverse();
	}
	void RotateFaces(FaceData& s){
		auto it = face_edge_start.begin();
		for (auto f: s){
			std::rotate(f->edges.begin(), f->edges.begin() + *it++, f->edges.end());
		}
	}
	void PermuteFaces(FaceData& s){
		ShpVector<Face> newf; newf.reserve(s.size());
		for (int i: face_order) newf.push_back(s[i]);
		std::swap(s, newf);
	}
public:
	RearrangeS(const FaceData& s, const Edge* ed0): input_surface(&s){
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

	void Apply(FaceData& s){
		ReverseEdges(s);
		RotateFaces(s);
		PermuteFaces(s);
	}
};

}//anonymous namespace

void hs::FaceRearrange(FaceData& s, const Edge* ed){
	RearrangeS r(s, ed);
	r.Apply(s);
}

// ================= Algos
bool hs::MatchTopology(const FaceData& a, const FaceData& b){
	if (a.size() != b.size()) return false;
	auto av = AllVertices(a);
	auto bv = AllVertices(b);
	if (av.size() != bv.size()) return false;

	//face->vertices connectivity
	std::vector<ShpVector<Vertex>> face_vert1, face_vert2;
	for (auto f: a) face_vert1.push_back(f->sorted_vertices());
	for (auto f: b) face_vert2.push_back(f->sorted_vertices());

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

EdgeData hs::ExtractBoundary(const FaceData& a, Vertex v){
	//leave only edges which have single face connection
	ShpVector<Edge> bedges;
	for (auto& fe: Connectivity::EdgeFace(a)){
		if (fe.size() == 1) bedges.push_back(fe.e);
	}
	return Connect(bedges, v);
}
std::array<vector<EdgeData>, 2> hs::ExtractAllBoundaries(const FaceData& a, Vect3 right_normal){
	std::array<vector<EdgeData>, 2> ret;
	//leave only boundary edges
	EdgeData bedges;
	vector<bool> reversed;
	for (auto& fe: Connectivity::EdgeFaceExtended(a)){
		if (fe.size() > 1) continue;
		bedges.push_back(fe.e);
		reversed.push_back(!fe.posdir[0]);
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
			} while (istart != icur && icur>=0 && used[icur]==false);

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

vector<FaceData> hs::ExtractSmooth(const FaceData& s, double angle){
	aa::enumerate_ids_pvec(s);
	vector<Vect3> normals(s.size());
	for (int i=0; i<s.size(); ++i) normals[i] = s[i]->left_normal();

	std::list<shared_ptr<Face>> unused(s.begin(), s.end());
	vector<FaceData> ret;

	auto use_face = [&](std::list<shared_ptr<Face>>::iterator it, FaceData* srf){
		(*srf).push_back(*it);
		return unused.erase(it);
	};

	double badcos = cos(angle*M_PI/180.);
	while (unused.size()>0){
		ret.push_back(FaceData());
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

	vector<FaceData> ret1;
	for (auto& r: ret){
		vector<FaceData> ss = Surface::AllSubSurfaces(r);
		ret1.resize(ret1.size()+ss.size());
		for (int i=0; i<ss.size(); ++i){
			std::swap(ret1[ret1.size()-ss.size()+i], ss[i]);
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
