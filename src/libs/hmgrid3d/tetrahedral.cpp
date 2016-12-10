#include "tetrahedral.hpp"
#include "Gmsh.h"
#include "GModel.h"
#include "MQuadrangle.h"
#include "MTriangle.h"
#include "MPoint.h"
#include "MLine.h"
#include "nan_handler.h"
#include "hmtimer.hpp"
using namespace HMGrid3D::Mesher;

namespace{

struct TransitionalFace{
	TransitionalFace(){}
	TransitionalFace(int lc, const std::vector<MVertex*>& vs){
		points = vs;
		for (auto p: points) ptsnums.insert(p->getIndex());
		left_cell = lc;
		right_cell = -1;
	}
	std::set<int> ptsnums;
	mutable int left_cell, right_cell;
	std::vector<MVertex*> points;
	mutable shared_ptr<HMGrid3D::Face> gface;
};

bool operator<(const TransitionalFace& f1, const TransitionalFace& f2){
	//comapare using only vertex indices
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
		if (p1->getIndex()>p2->getIndex()) std::swap(p1, p2);
	}
	MVertex *p1, *p2;
	int ind;
};
bool operator<(const TransitionalEdge& e1, const TransitionalEdge& e2){
	if (e1.p1->getIndex() < e2.p1->getIndex()) return true;
	if (e1.p1->getIndex() > e2.p1->getIndex()) return false;
	return (e1.p2->getIndex() < e2.p2->getIndex());
}

std::set<TransitionalFace> GridFromModel(GModel& m, HMGrid3D::GridData& ret){
	ret.clear();
	m.indexMeshVertices(true);
	ret.vvert.resize(m.getNumMeshVertices());
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
				int index = v->getIndex()-1;
				assert(index<ret.vvert.size());
				if (!ret.vvert[index]) ret.vvert[index].reset(new HMGrid3D::Vertex(v->x(), v->y(), v->z()));
			}
			//faces to transitional status:: Takes 70% of exec time
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
				throw std::runtime_error("unsupported gmsh element type");
			}
			aa::add_shared(ret.vcells, HMGrid3D::Cell());
		}
	}
	//faces && edges: Takes 40% of exec time;
	int numed = 0;
	for (auto& f: transface){
		auto nf = aa::add_shared(ret.vfaces, HMGrid3D::Face());
		f.gface=ret.vfaces.back();
		for (int i=0; i<f.ptsnums.size(); ++i){
			int inext = (i+1)%f.ptsnums.size();
			auto m1 = f.points[i], m2 = f.points[inext];
			auto te = transedge.emplace(m1, m2, numed);
			if (te.second) {
				auto ind1 = m1->getIndex()-1;
				auto ind2 = m2->getIndex()-1;
				aa::add_shared(ret.vedges, HMGrid3D::Edge(ret.vvert[ind1], ret.vvert[ind2]));
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
	return transface;
}

GFace* fill_model_with_2d(GModel& m, const HMGrid3D::Surface& surf, vector<vector<GEdge*>>& gedges, vector<MVertex*>& mvheap){
	//build gmsh faces
	GFace* g_face = m.addPlanarFace(gedges);

	auto place_mvertex = [&](HMGrid3D::Vertex& v)->MVertex*{
		if (mvheap[v.id] == 0){
			mvheap[v.id] = new MVertex(v.x, v.y, v.z, g_face);
			g_face->mesh_vertices.push_back(mvheap[v.id]);
		}
		return mvheap[v.id];
	};
	//add mesh
	for (int i=0; i<surf.faces.size(); ++i){
		auto sv = surf.faces[i]->sorted_vertices();
		vector<MVertex*> mv(sv.size());

		for (int j=0; j<sv.size(); ++j) mv[j] = place_mvertex(*sv[j]);

		if (sv.size()==4){
			auto nq = new MQuadrangle(mv[0], mv[1], mv[2], mv[3]);
			g_face->addQuadrangle(nq);
		} else if (sv.size()==3){
			auto nq = new MTriangle(mv[0], mv[1], mv[2]);
			g_face->addTriangle(nq);
		} else {
			throw std::runtime_error("Only triangle and quadrangle surface faces are allowed");
		};
	}

	return g_face;
}
vector<vector<GEdge*>> fill_model_with_1d(
		GModel& m,
		const vector<HMGrid3D::EdgeData>& edata,
		vector<GEdge*>& eheap, vector<GVertex*>& vheap,
		vector<MVertex*>& mvheap){
	auto place_vertex = [&](HMGrid3D::Vertex& v, double dist)->GVertex*{
		if (vheap[v.id] == 0){
			vheap[v.id] = m.addVertex(v.x, v.y, v.z, dist);
			mvheap[v.id] = new MVertex(v.x, v.y, v.z, vheap[v.id]);
			vheap[v.id]->mesh_vertices.push_back(mvheap[v.id]);
			vheap[v.id]->points.push_back(new MPoint(mvheap[v.id]));
		}
		auto ms = vheap[v.id]->prescribedMeshSizeAtVertex();
		if (ms>dist) vheap[v.id]->setPrescribedMeshSizeAtVertex(dist);
		return vheap[v.id];
	};
	auto place_edge = [&](HMGrid3D::Edge& ecur, HMGrid3D::Edge& eprev)->GEdge*{
		if (eheap[ecur.id] != 0) return eheap[ecur.id];
		auto v1 = ecur.first().get();
		auto v2 = ecur.last().get();
		if (v1 != eprev.first().get() && v1 != eprev.last().get()){
			std::swap(v1, v2);
		}
		double dist = HMGrid3D::Vertex::dist(*v1, *v2);
		auto gv1 = place_vertex(*v1, dist);
		auto gv2 = place_vertex(*v2, dist);
		eheap[ecur.id] = m.addLine(gv1, gv2);
		eheap[ecur.id]->addLine(new MLine(mvheap[v1->id], mvheap[v2->id]));
		return eheap[ecur.id];
	};
	vector<vector<GEdge*>> ret;
	for (auto& ed: edata){
		ret.emplace_back();
		auto& r = ret.back();
		HMGrid3D::Edge* prev = ed.back().get();
		for (auto& e: ed){
			r.push_back(place_edge(*e, *prev));
			prev = e.get();
		}
	}
	return ret;
}

void fill_model_with_3d(GModel& m, const vector<vector<GFace*>>& fc){
	//Mesh3D
	auto volume = m.addVolume(fc);
	//m.writeGEO("gmsh_geo.geo");
	//m.writeMSH("gmsh_geo.msh");
	NanSignalHandler::StopCheck();
	m.mesh(3);
	NanSignalHandler::StartCheck();
	//m.writeVTK("gmsh_geo.vtk");
	//m.writeMSH("gmsh_geo.msh");
}

void restore_btypes(GFace* gf, std::set<TransitionalFace>& fcset, const HMGrid3D::Surface& surf){
	auto analyze_melement = [&fcset](MElement* el, int bt){
		TransitionalFace fndface;
		for (int i=0; i<el->getNumVertices(); ++i){
			fndface.ptsnums.insert(
				el->getVertex(i)->getIndex()
			);
		}
		auto fnd = fcset.find(fndface);
		if (fnd != fcset.end()){
			fnd->gface->boundary_type = bt;
			fcset.erase(fnd);
		} else {assert(false);}
	};

	//keeping the same order at which 2d grid for this surface was created
	//add mesh
	int itris=-1, iquads=-1;
	for (int i=0; i<surf.faces.size(); ++i){
		auto& f = surf.faces[i];
		if (f->n_edges() == 3) ++itris;
		else if (f->n_edges() == 4) ++iquads;
		else assert(false);

		int btype = f->boundary_type;
		if (btype==0) continue;
		
		auto sv = surf.faces[i]->sorted_vertices();
		vector<MVertex*> mv(sv.size());

		MElement* nq;
		if (sv.size()==4) nq = gf->quadrangles[iquads];
		else if (sv.size()==3) nq = gf->triangles[itris];
		else throw std::runtime_error("Only triangle and quadrangle surface faces are allowed");

		for (int j=0; j<sv.size(); ++j) mv[j] = nq->getVertex(j);
		analyze_melement(nq, btype);
	}
}

HMGrid3D::GridData gmsh_fill(const HMGrid3D::SurfaceTree& tree, const HMGrid3D::Surface& cond,
		const HMGrid3D::VertexData& pcond, const vector<double>& psizes,
		HMCallback::Caller2& cb){
	assert(
		//first node is root, all others are its children
		[&](){
			if (tree.nodes[0]->level != 0) return false;
			for (int i=1; i<tree.nodes.size(); ++i)
				if (tree.nodes[i]->level != 1) return false;
			return true;
		}()
	);
	GModel m;
	m.setFactory("Gmsh");
	GmshSetOption("General", "Verbosity", 0.0);
	//GmshSetOption("Mesh", "OptimizeNetgen", 1.0);
	
	//decomposition
	cb.silent_step_after(10, "Surfaces decomposition");
	HMGrid3D::EdgeData ae=tree.alledges();
	HMGrid3D::VertexData av=tree.allvertices();
	vector<vector<HMGrid3D::Surface>> decomposed_surfs;
	for (auto n: tree.nodes){
		decomposed_surfs.push_back(HMGrid3D::Surface::ExtractSmooth(*n, 30));
	}
	vector<vector<vector<HMGrid3D::EdgeData>>> decomposed_edges;
	for (auto& ds: decomposed_surfs){
		vector<vector<HMGrid3D::EdgeData>> de;
		for (auto& s: ds){
			vector<HMGrid3D::EdgeData> ds;
			Vect3 right_normal = s.faces[0]->left_normal()*(-1);
			auto ex = HMGrid3D::Surface::ExtractAllBoundaries(s, right_normal);
			assert(ex[0].size() == 1);
			ds.push_back(ex[0][0]);
			for (int i=0; i<ex[1].size(); ++i) ds.push_back(ex[1][i]);
			de.push_back(ds);
		}
		decomposed_edges.push_back(de);
	}


	//mesh1d
	cb.step_after(5, "Fill 1D mesh");
	enumerate_ids_pvec(ae);
	enumerate_ids_pvec(av);
	vector<GEdge*> g_edges_heap(ae.size(), 0);
	vector<GVertex*> g_vertex_heap(av.size(), 0);
	vector<MVertex*> m_vertex_heap(av.size(), 0);
	vector<vector<vector<GEdge*>>> g_edges;
	for(int id=0; id<decomposed_surfs.size(); ++id)
	for(int is=0; is<decomposed_surfs[id].size(); ++is){
		auto& d = decomposed_edges[id][is];
		g_edges.push_back(fill_model_with_1d(m, d, g_edges_heap, g_vertex_heap, m_vertex_heap));
	}

	//mesh2d
	cb.step_after(10, "Fill 2D mesh");
	vector<vector<GFace*>> g_faces;
	int k=0;
	for (auto& ds: decomposed_surfs){
		vector<GFace*> g_faces_1;
		for (auto& s: ds){
			GFace* gf = fill_model_with_2d(m, s, g_edges[k++], m_vertex_heap);
			g_faces_1.push_back(gf);
		}
		g_faces.push_back(g_faces_1);
	}

	//mesh3d
	cb.step_after(35, "Build 3D mesh");
	fill_model_with_3d(m, g_faces);

	HMGrid3D::GridData ret;
	cb.step_after(30, "Assemble mesh");
	auto tfaces = GridFromModel(m, ret);

	//restore boundary types from tree
	cb.step_after(10, "Restore boundary types");
	auto fit = tfaces.begin();
	while (fit!=tfaces.end()){
		if (fit->left_cell>=0 && fit->right_cell>=0){
			fit = tfaces.erase(fit);
		} else ++fit;
	}
	ret.enumerate_all();
	int k1=0;
	for (auto& ds: decomposed_surfs){
		vector<GFace*>& g_faces_1 = g_faces[k1++];
		int k2=0;
		for (auto& s: ds){
			GFace* gf = g_faces_1[k2++];
			restore_btypes(gf, tfaces, s);
		}
	}

	cb.fin();
	return ret;
}
};

HMCallback::FunctionWithCallback<TUnstructuredTetrahedral> HMGrid3D::Mesher::UnstructuredTetrahedral;

HMGrid3D::GridData TUnstructuredTetrahedral::_run(const HMGrid3D::Surface& source,
		const HMGrid3D::Surface& sinner,
		const VertexData& pinner, const vector<double>& psizes){
	callback->step_after(20, "Surfaces nesting");
	//main tree
	HMGrid3D::SurfaceTree stree = HMGrid3D::SurfaceTree::Assemble(source);
	if (stree.roots().size() == 0) throw std::runtime_error("No closed surfaces found");
	//cropped trees
	vector<HMGrid3D::SurfaceTree> trees = stree.crop_level1();
	//conditional surfaces
	HMGrid3D::Surface cond;
	std::copy(sinner.faces.begin(), sinner.faces.end(), std::back_inserter(cond.faces));
	for (auto nd: stree.nodes) if (nd->isopen()){
		std::copy(nd->faces.begin(), nd->faces.end(), std::back_inserter(cond.faces));
	}
	//temporary revert faces to match surfaces
	ShpVector<HMGrid3D::SurfTReverterBase> revs;
	revs.emplace_back(new HMGrid3D::SurfTReverter(cond));
	for (auto tree: trees) revs.emplace_back(new HMGrid3D::SurfTreeTReverter(tree));

	HMGrid3D::GridData ret;
	for (int i=0; i<trees.size(); ++i){
		auto cb = callback->subrange(75./trees.size(), 100.);
		if (i == 0) ret = gmsh_fill(trees[i], cond, pinner, psizes, *cb);
		else {
			HMGrid3D::GridData sg = gmsh_fill(trees[i], cond, pinner, psizes, *cb);
			std::copy(sg.vvert.begin(), sg.vvert.end(), std::back_inserter(ret.vvert));
			std::copy(sg.vedges.begin(), sg.vedges.end(), std::back_inserter(ret.vedges));
			std::copy(sg.vfaces.begin(), sg.vfaces.end(), std::back_inserter(ret.vfaces));
			std::copy(sg.vcells.begin(), sg.vcells.end(), std::back_inserter(ret.vcells));
		}
	}
	callback->step_after(5, "Finalizing");
	return ret;
}

HMGrid3D::GridData TUnstructuredTetrahedral::_run(const HMGrid3D::Surface& source){
	return _run(source, HMGrid3D::Surface(), {}, {});
}
HMGrid3D::GridData TUnstructuredTetrahedral::_run(const HMGrid3D::Surface& source,
		const HMGrid3D::Surface& sinner){
	return _run(source, sinner, {}, {});
}
HMGrid3D::GridData TUnstructuredTetrahedral::_run(const HMGrid3D::Surface& source,
		const VertexData& pinner, const vector<double>& psizes){
	return _run(source, HMGrid3D::Surface(), pinner, psizes);
}
