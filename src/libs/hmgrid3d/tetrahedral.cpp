#include "tetrahedral.hpp"
#include "Gmsh.h"
#include "GModel.h"
#include "MQuadrangle.h"
#include "MTriangle.h"
#include "MElementCut.h"
#include "MPoint.h"
#include "MLine.h"
#include "nan_handler.h"
#include "hmtimer.hpp"
#include "tetramesh_preproc.hpp"
#include "debug3d.hpp"
#include "treverter3d.hpp"
using namespace HM3D::Mesher;
using namespace HM3D;

namespace{

struct TransitionalFace{
	TransitionalFace(){}
	TransitionalFace(TransitionalFace&&) = default;
	TransitionalFace(const TransitionalFace&) = default;
	TransitionalFace& operator=(const TransitionalFace&) = default;

	TransitionalFace(int lc, const std::vector<MVertex*>& vs):points(vs), ptsnums(vs.size()){
		for (int i=0; i<points.size(); ++i) ptsnums[i] = points[i]->getIndex();
		std::sort(ptsnums.begin(), ptsnums.end());
		left_cell = lc;
		right_cell = -1;
	}
	TransitionalFace(int lc, MVertex* v1, MVertex* v2, MVertex* v3):points(3), ptsnums(3){
		points[0] = v1; ptsnums[0] = v1->getIndex();
		points[1] = v2; ptsnums[1] = v2->getIndex();
		points[2] = v3; ptsnums[2] = v3->getIndex();
		std::sort(ptsnums.begin(), ptsnums.end());
		left_cell = lc;
		right_cell = -1;
	}
	std::vector<int> ptsnums;
	mutable int left_cell, right_cell;
	std::vector<MVertex*> points;
	mutable shared_ptr<HM3D::Face> gface;
};

bool operator<(const TransitionalFace& f1, const TransitionalFace& f2){
	//compare using only vertex indices
	if (f1.ptsnums[0] < f2.ptsnums[0]) return true;
	else if (f1.ptsnums[0] > f2.ptsnums[0]) return false;
	else if (f1.ptsnums[1] < f2.ptsnums[1]) return true;
	else if (f1.ptsnums[1] > f2.ptsnums[1]) return false;
	return (f1.ptsnums[2] < f2.ptsnums[2]);
}
bool operator!=(const TransitionalFace& f1, const TransitionalFace& f2){
	//compare using only vertex indices
	if (f1.ptsnums.size() != f2.ptsnums.size()) return true;
	auto it1 = f1.ptsnums.begin(), it2 = f2.ptsnums.begin();
	while (it1 != f1.ptsnums.end()){
		if (*it1 != *it2) return true;
		++it1; ++it2;
	};
	return false;
}
struct TransitionalEdge{
	TransitionalEdge(MVertex* a, MVertex* b, int num):
			p1(a), p2(b), ind(num),
			i1(p1->getIndex(), p2->getIndex()){
		if (i1.first>i1.second) std::swap(i1.first, i1.second);
	}
	MVertex *p1, *p2;
	int ind;
	std::pair<int, int> i1;
};
bool operator<(const TransitionalEdge& e1, const TransitionalEdge& e2){
	return e1.i1 < e2.i1;
}

void GridFromModel(GModel& m, HM3D::GridData& ret){
	ret.clear();
	m.indexMeshVertices(true);
	ret.vvert.resize(m.getNumMeshVertices());
	std::set<TransitionalFace> transface;
	std::set<TransitionalEdge> transedge;

	auto add_face = [&](const vector<MVertex*>& elvert, int en){
		auto tf = transface.emplace(en, elvert);
		if (tf.second == false) tf.first->right_cell = en;
	};
	auto add_face3 = [&](MVertex* v1, MVertex* v2, MVertex* v3, int en){
		auto tf = transface.emplace(en, v1, v2, v3);
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
				if (!ret.vvert[index]) ret.vvert[index].reset(new HM3D::Vertex(v->x(), v->y(), v->z()));
			}
			//faces to transitional status:: Takes 50% of function exec time
			if (e->getType() == TYPE_TET){
				add_face3(elvert[0], elvert[2], elvert[1], en);
				add_face3(elvert[0], elvert[3], elvert[2], en);
				add_face3(elvert[1], elvert[2], elvert[3], en);
				add_face3(elvert[0], elvert[1], elvert[3], en);
			} else if (e->getType() == TYPE_PYR){
				add_face3(elvert[0], elvert[1], elvert[4], en);
				add_face3(elvert[1], elvert[2], elvert[4], en);
				add_face3(elvert[2], elvert[3], elvert[4], en);
				add_face3(elvert[0], elvert[4], elvert[3], en);
				add_face({elvert[0], elvert[3], elvert[2], elvert[1]}, en);
			} else {
				throw std::runtime_error("unsupported gmsh element type");
			}
			aa::add_shared(ret.vcells, HM3D::Cell());
		}
	}
	while (ret.vvert.size() > 0 && ret.vvert.back() == nullptr)
		ret.vvert.resize(ret.vvert.size()-1);
	//check if grid was ok => all vertices were read
	if (ret.vvert.size() == 0 || std::any_of(ret.vvert.begin(), ret.vvert.end(),
			[](shared_ptr<HM3D::Vertex>& v){return v==nullptr;})){
		throw std::runtime_error("Tetrahedral meshing failed");
	}

	//faces && edges: Takes 50% of function exec time;
	int numed = 0;
	for (auto& f: transface){
		auto nf = aa::add_shared(ret.vfaces, HM3D::Face());
		f.gface=ret.vfaces.back();
		for (int i=0; i<f.ptsnums.size(); ++i){
			int inext = (i+1)%f.ptsnums.size();
			auto m1 = f.points[i], m2 = f.points[inext];
			auto te = transedge.emplace(m1, m2, numed);
			if (te.second) {
				auto ind1 = m1->getIndex()-1;
				auto ind2 = m2->getIndex()-1;
				ret.vedges.emplace_back(new HM3D::Edge(ret.vvert[ind1], ret.vvert[ind2]));
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
}

GFace* fill_model_with_2d(GModel& m, const HM3D::FaceData& surf, vector<vector<GEdge*>>& gedges, vector<MVertex*>& mvheap){
	//build gmsh faces
	GFace* g_face = m.addPlanarFace(gedges);

	auto place_mvertex = [&](HM3D::Vertex& v)->MVertex*{
		if (mvheap[v.id] == 0){
			mvheap[v.id] = new MVertex(v.x, v.y, v.z, g_face);
			g_face->mesh_vertices.push_back(mvheap[v.id]);
		}
		return mvheap[v.id];
	};
	//add mesh
	for (int i=0; i<surf.size(); ++i){
		auto sv = surf[i]->sorted_vertices();
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
		const vector<HM3D::EdgeData>& edata,
		vector<GEdge*>& eheap, vector<GVertex*>& vheap,
		vector<MVertex*>& mvheap){
	auto place_vertex = [&](HM3D::Vertex& v, double dist)->GVertex*{
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
	auto place_edge = [&](HM3D::Edge& ecur, HM3D::Edge& eprev)->GEdge*{
		if (eheap[ecur.id] != 0) return eheap[ecur.id];
		auto v1 = ecur.first().get();
		auto v2 = ecur.last().get();
		if (v1 != eprev.first().get() && v1 != eprev.last().get()){
			std::swap(v1, v2);
		}
		double dist = HM3D::Vertex::dist(*v1, *v2);
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
		HM3D::Edge* prev = ed.back().get();
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
	auto bb = m.bounds();
	GmshSetBoundingBox(bb.min()[0], bb.max()[0], bb.min()[1], bb.max()[1], bb.min()[2], bb.max()[2]);
	NanSignalHandler::StopCheck();
	m.mesh(3);
	NanSignalHandler::StartCheck();
	//m.writeVTK("gmsh_geo.vtk");
	//m.writeMSH("gmsh_geo.msh");
}

void equal_vertices(const HM3D::GridData& grid, GFace* gf, const FaceData& surf,
		HM3D::VertexData& gfvert, HM3D::VertexData& svert){
	//here we rely on fact that GFace mesh vertices are indexed according to grid.vvert order.
	//this is guaranteed by earlier GridFromModel procedure call
	HM3D::VertexData s = AllVertices(surf);
	HM3D::VertexData g2;
	auto addvert = [&g2, &grid](MElement* el){
		for (int j=0; j<el->getNumVertices(); ++j){
			g2.push_back(grid.vvert[el->getVertex(j)->getIndex() - 1]);
		}
	};
	for (auto nq: gf->triangles) addvert(nq);
	for (auto nq: gf->quadrangles) addvert(nq);
	for (auto nq: gf->polygons) addvert(nq);

	HM3D::VertexData g;
	aa::constant_ids_pvec(g2, -1);
	for (auto it: g2) if (it->id == -1){
		g.push_back(it);
		it->id = 0;
	}
	assert(s.size() == g.size());
	auto vertsort = [](shared_ptr<HM3D::Vertex> a, shared_ptr<HM3D::Vertex> b)->bool{
		if (ISLOWER(a->x, b->x)) return true;
		else if (ISLOWER(b->x, a->x)) return false;
		else if (ISLOWER(a->y, b->y)) return true;
		else if (ISLOWER(b->y, a->y)) return false;
		else return ISLOWER(a->z, b->z);
	};
	std::sort(s.begin(), s.end(), vertsort);
	std::sort(g.begin(), g.end(), vertsort);
	gfvert.resize(gfvert.size()+g.size());
	svert.resize(svert.size()+s.size());
	std::copy(g.begin(), g.end(), gfvert.end() - g.size());
	std::copy(s.begin(), s.end(), svert.end() - s.size());
}

HM3D::GridData gmsh_fill(const Surface::Tree& tree, const FaceData& cond,
		const HM3D::VertexData& pcond, const vector<double>& psizes,
		HMCallback::Caller2& cb){
	GModel m;
	m.setFactory("Gmsh");
	GmshSetOption("General", "Verbosity", 0.0);
	GmshSetOption("Mesh", "Optimize", 1.0);
	//GmshSetOption("Mesh", "OptimizeNetgen", 1.0);
	
	//decomposition
	cb.step_after(10, "Boundary preprocessing");
	HM3D::Mesher::SurfacePreprocess presurf(tree, 30);
	//if whole area is meshed with pyramids
	if (presurf.decomposed_surfs.size() == 0){
		cb.fin();
		return presurf.bnd_grid;
	}
	
	//mesh1d
	cb.step_after(5, "Fill 1D mesh");
	aa::enumerate_ids_pvec(presurf.ae);
	aa::enumerate_ids_pvec(presurf.av);
	vector<GEdge*> g_edges_heap(presurf.ae.size(), 0);
	vector<GVertex*> g_vertex_heap(presurf.av.size(), 0);
	vector<MVertex*> m_vertex_heap(presurf.av.size(), 0);
	vector<vector<vector<GEdge*>>> g_edges;
	for(int id=0; id<presurf.decomposed_surfs.size(); ++id)
	for(int is=0; is<presurf.decomposed_surfs[id].size(); ++is){
		auto& d = presurf.decomposed_edges[id][is];
		g_edges.push_back(fill_model_with_1d(m, d, g_edges_heap, g_vertex_heap, m_vertex_heap));
	}

	//mesh2d
	cb.step_after(10, "Fill 2D mesh");
	vector<vector<GFace*>> g_faces;
	int k=0;
	for (auto& ds: presurf.decomposed_surfs){
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

	HM3D::GridData ret;
	cb.step_after(30, "Assemble mesh");
	GridFromModel(m, ret);

	//restore boundary types from tree
	cb.step_after(10, "Restore boundary");
	HM3D::VertexData r1, s1;
	int k1=0;
	for (auto& ds: presurf.decomposed_surfs){
		vector<GFace*>& g_faces_1 = g_faces[k1++];
		int k2=0;
		for (auto& s: ds){
			GFace* gf = g_faces_1[k2++];
			equal_vertices(ret, gf, s, r1, s1);
		}
	}

	HM3D::VertexData ret_duplicates;
	HM3D::VertexData surfs_duplicates;
	aa::constant_ids_pvec(r1, -1);
	for (int i=0; i<r1.size(); ++i) if (r1[i]->id == -1){
		r1[i]->id = 0;
		ret_duplicates.push_back(r1[i]);
		surfs_duplicates.push_back(s1[i]);
	}
	presurf.Restore(ret, ret_duplicates, surfs_duplicates);

	cb.fin();
	return ret;
}
};

HMCallback::FunctionWithCallback<TUnstructuredTetrahedral> HM3D::Mesher::UnstructuredTetrahedral;

HM3D::GridData TUnstructuredTetrahedral::_run(const FaceData& source,
		const FaceData& sinner,
		const VertexData& pinner, const vector<double>& psizes){
	callback->step_after(20, "Surfaces nesting");
	//main tree
	Surface::Tree stree = Surface::Tree::Assemble(source);
	if (stree.roots().size() == 0) throw std::runtime_error("No closed surfaces found");
	//cropped trees
	vector<Surface::Tree> trees = stree.crop_level1();
	//conditional surfaces
	FaceData cond;
	if (sinner.size() > 0) _THROW_NOT_IMP_;
	/* TODO
	std::copy(sinner.faces.begin(), sinner.faces.end(), std::back_inserter(cond.faces));
	for (auto nd: stree.nodes) if (nd->isopen()){
		std::copy(nd->faces.begin(), nd->faces.end(), std::back_inserter(cond.faces));
	}
	*/
	//temporary revert faces to match surfaces
	ShpVector<Surface::RevertTree> revs;
	//revs.emplace_back(new HM3D::SurfTReverter(cond));
	for (auto tree: trees) revs.emplace_back(new Surface::RevertTree(tree));

	HM3D::GridData ret;
	for (int i=0; i<trees.size(); ++i){
		auto cb = callback->subrange(75./trees.size(), 100.);
		if (i == 0) ret = gmsh_fill(trees[i], cond, pinner, psizes, *cb);
		else {
			HM3D::GridData sg = gmsh_fill(trees[i], cond, pinner, psizes, *cb);
			std::copy(sg.vvert.begin(), sg.vvert.end(), std::back_inserter(ret.vvert));
			std::copy(sg.vedges.begin(), sg.vedges.end(), std::back_inserter(ret.vedges));
			std::copy(sg.vfaces.begin(), sg.vfaces.end(), std::back_inserter(ret.vfaces));
			std::copy(sg.vcells.begin(), sg.vcells.end(), std::back_inserter(ret.vcells));
		}
	}
	
	callback->step_after(5, "Finalizing");

	//check volumes: this check should be abandoned
	//when algorithm becomes more reliable
	double v1 = Surface::Volume(source);
	double v2 = 0.;
	for (auto& c: ret.vcells){
		double vc = c->volume();
		assert(vc > 0);
		v2 += vc;
	}
	if (!ISEQ(v1, v2)) throw std::runtime_error(
		"Volumes check failed");

	
	return ret;
}

HM3D::GridData TUnstructuredTetrahedral::_run(const FaceData& source){
	return _run(source, FaceData(), {}, {});
}
HM3D::GridData TUnstructuredTetrahedral::_run(const FaceData& source,
		const FaceData& sinner){
	return _run(source, sinner, {}, {});
}
HM3D::GridData TUnstructuredTetrahedral::_run(const FaceData& source,
		const VertexData& pinner, const vector<double>& psizes){
	return _run(source, FaceData(), pinner, psizes);
}
