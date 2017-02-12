#include "assemble3d.hpp"
#include "contabs3d.hpp"
#include "finder3d.hpp"
using namespace HM3D;

namespace hs = HM3D::Surface::Assembler;
namespace hc = HM3D::Contour::Assembler;

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

FaceData hs::SubSurface(const FaceData& s, const Vertex* v){
	//subdivide
	vector<FaceData> subd = HM3D::SplitData(s);
	//find part which contains v
	for (auto& fvec: subd){
		auto av = AllVertices(fvec);
		if (aa::shp_find(av.begin(), av.end(), v) != av.end()) return fvec;
	}
	return FaceData();
}


EdgeData hc::ExtractBoundary(const FaceData& a, Vertex v){
	//leave only edges which have single face connection
	ShpVector<Edge> bedges;
	for (auto& fe: Connectivity::EdgeFace(a)){
		if (fe.size() == 1) bedges.push_back(fe.e);
	}
	return Connect(bedges, v);
}
std::array<vector<EdgeData>, 2> hc::ExtractAllBoundaries(const FaceData& a, Vect3 right_normal){
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
		vector<FaceData> ss = HM3D::SplitData(r);
		ret1.resize(ret1.size()+ss.size());
		for (int i=0; i<ss.size(); ++i){
			std::swap(ret1[ret1.size()-ss.size()+i], ss[i]);
		}
	}
	
	return ret1;
}

EdgeData hc::Connect(const EdgeData& data, Vertex app_v){
	//1. find closest vertex to v
	ShpVector<Vertex> vall;
	for (auto d: data) { vall.push_back(d->first()); vall.push_back(d->last()); }
	vall = aa::no_duplicates(vall);
	auto fres = Finder::ClosestPoint(vall, app_v);
	const shared_ptr<Vertex> v = vall[std::get<0>(fres)];

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
