#include "tetrahedral.hpp"
#include "GModel.h"
#include "MQuadrangle.h"
#include "MTriangle.h"
#include "MPoint.h"
#include "MLine.h"
#include "nan_handler.h"
using namespace HMGrid3D;

template<class A>
ShpVector<A> enumerate_unique(const vector<ShpVector<A>>& input){
	for (auto a: input)
	for (auto b: a) b->id = -1;

	ShpVector<A> ret;
	for (auto a: input)
	for (auto b: a) if (b->id == -1){
		b->id = ret.size();
		ret.push_back(b);
	}

	return ret;
}

namespace{
struct TransitionalFace{
	TransitionalFace(int lc, const std::vector<MVertex*>& vs){
		points = vs;
		for (auto p: points) ptsnums.insert(p->getNum());
		left_cell = lc;
		right_cell = -1;
	}
	std::set<int> ptsnums;
	mutable int left_cell, right_cell;
	std::vector<MVertex*> points;
};

bool operator<(const TransitionalFace& f1, const TransitionalFace& f2){
	if (f1.ptsnums.size() < f2.ptsnums.size()) return true;
	if (f1.ptsnums.size() > f2.ptsnums.size()) return false;
	auto it1 = f1.ptsnums.begin(), it2 = f2.ptsnums.begin();
	while (it1 != f1.ptsnums.end()){
		if (*it1 < *it2) return true;
		if (*it1 > *it2) return false;
		++it1; ++it2;
	};
	return false;
}

struct TransitionalEdge{
	TransitionalEdge(MVertex* a, MVertex* b, int num): p1(a), p2(b), ind(num){
		if (p1->getNum()>p2->getNum()) std::swap(p1, p2);
	}
	MVertex *p1, *p2;
	int ind;
};
bool operator<(const TransitionalEdge& e1, const TransitionalEdge& e2){
	if (e1.p1->getNum() < e2.p1->getNum()) return true;
	if (e1.p1->getNum() > e2.p1->getNum()) return false;
	return (e1.p2->getNum() < e2.p2->getNum());
}

};

SGrid GridFromModel(GModel& m){
	SGrid ret;

	std::map<int, shared_ptr<Vertex>> vertmap;
	std::set<TransitionalFace> transface;
	std::set<TransitionalEdge> transedge;
	auto add_face = [&](const vector<MVertex*>& elvert, int en){
		auto tf = transface.emplace(en, elvert);
		if (tf.second == false) tf.first->right_cell = en;
	};
	for (auto it = m.firstRegion(); it!=m.lastRegion(); ++it){
		for (int en=0; en<(*it)->getNumMeshElements(); ++en){
			auto e = (*it)->getMeshElement(en);
			//vertices
			std::vector<MVertex*> elvert;
			e->getVertices(elvert);
			for (int i=0; i<e->getNumVertices(); ++i){
				auto v = elvert[i];
				auto empres = vertmap.emplace(v->getNum(), shared_ptr<Vertex>());
				if (empres.second){
					empres.first->second.reset(new Vertex(v->x(), v->y(), v->z()));
					ret.vvert.push_back(empres.first->second);
				}
			}
			
			//faces to transitional status
			if (e->getType() == TYPE_TET){
				add_face({elvert[0], elvert[2], elvert[1]}, en);
				add_face({elvert[0], elvert[3], elvert[2]}, en);
				add_face({elvert[1], elvert[2], elvert[3]}, en);
				add_face({elvert[0], elvert[1], elvert[3]}, en);
			} else if (e->getType() == TYPE_PYR){
				add_face({elvert[0], elvert[1], elvert[4]}, en);
				add_face({elvert[1], elvert[2], elvert[4]}, en);
				add_face({elvert[2], elvert[3], elvert[4]}, en);
				add_face({elvert[0], elvert[4], elvert[3]}, en);
				add_face({elvert[0], elvert[3], elvert[2], elvert[1]}, en);
			} else {
				throw std::runtime_error("unknown element type");
			}
			aa::add_shared(ret.vcells, HMGrid3D::Cell());
		}
	}
	//faces && edges
	int numed = 0;
	for (auto& f: transface){
		auto nf = aa::add_shared(ret.vfaces, HMGrid3D::Face());
		for (int i=0; i<f.ptsnums.size(); ++i){
			int inext = (i+1)%f.ptsnums.size();
			auto m1 = f.points[i], m2 = f.points[inext];
			auto te = transedge.emplace(m1, m2, numed);
			if (te.second) {
				aa::add_shared(ret.vedges, HMGrid3D::Edge(vertmap[m1->getNum()], vertmap[m2->getNum()]));
				++numed;
			}
			nf->edges.push_back(ret.vedges[te.first->ind]);
		}
		if (f.left_cell>=0) nf->left = ret.vcells[f.left_cell];
		if (f.right_cell>=0) nf->right = ret.vcells[f.right_cell];
	}
	//cell ->face connectivity
	for (auto& f: ret.vfaces){
		if (f->has_left_cell()) f->left.lock()->faces.push_back(f);
		if (f->has_right_cell()) f->right.lock()->faces.push_back(f);
	}
	ret.actualize_serial_data();
	return ret;
}

SGrid Mesher::UnstructedTetrahedral(const ShpVector<HMGrid3D::Surface>& srfs){
	vector<ShpVector<Edge>> f_alledges;
	vector<ShpVector<Vertex>> f_allvertices;
	
	for (auto s: srfs){
		f_alledges.push_back(s->alledges());
		f_allvertices.push_back(s->allvertices());
	}
	
	ShpVector<Edge> alledges = enumerate_unique(f_alledges);
	ShpVector<Vertex> allvertices = enumerate_unique(f_allvertices);

	GModel m;
	m.setFactory("Gmsh");

	vector<GVertex*> g_vertex;
	for (auto v: allvertices) g_vertex.push_back(m.addVertex(v->x, v->y, v->z, 0.2));
	vector<GEdge*> g_edge;
	for (auto e: alledges){
		int p1 = e->first()->id;
		int p2 = e->last()->id;
		g_edge.push_back(m.addLine(g_vertex[p1], g_vertex[p2]));
	}

	vector<GFace*> g_face;
	ShpVector<Face> allfaces;
	for (auto s: srfs){
		for (auto f: s->faces){
			allfaces.push_back(f);
			vector<GEdge*> eds;
			for (auto e: f->edges){
				eds.push_back(g_edge[e->id]);
			}
			g_face.push_back(m.addPlanarFace({eds}));
		}
	}
	auto volume = m.addVolume({g_face});

	//Mesh0D
	vector<MVertex*> mvertex;
	for (int i=0; i<g_vertex.size(); ++i){
		GVertex *gv = g_vertex[i];
		gv->mesh_vertices.push_back(new MVertex(gv->x(), gv->y(), gv->z(), gv));
		mvertex.push_back(gv->mesh_vertices.back());
		gv->points.push_back(new MPoint(gv->mesh_vertices.back()));
	}
        //Mesh1D
	vector<MEdge*> medge;
	for (int i=0; i<g_edge.size(); ++i){
		int p1 = alledges[i]->first()->id;
		int p2 = alledges[i]->last()->id;
		g_edge[i]->addLine(new MLine(mvertex[p1], mvertex[p2]));
	}
	//Mesh2D
	vector<MQuadrangle*> mquads;
	for (int i=0; i<allfaces.size(); ++i){
		auto sv = allfaces[i]->sorted_vertices();
		auto nq = new MQuadrangle(
			mvertex[sv[0]->id],
			mvertex[sv[1]->id],
			mvertex[sv[2]->id],
			mvertex[sv[3]->id]);
		g_face[i]->addQuadrangle(nq);
	}

	//Mesh3D
	//m.writeGEO("gmsh_geo.geo");
	NanSignalHandler::StopCheck();
	m.mesh(3);
	NanSignalHandler::StartCheck();
	m.writeVTK("gmsh_geo.vtk");
	m.writeMSH("gmsh_geo.msh");

	return GridFromModel(m);
}
