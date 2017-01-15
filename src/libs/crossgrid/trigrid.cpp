#include "trigrid.hpp"
#include "algos.hpp"
#include "cont_partition.hpp"
#include "contabs2d.hpp"
#include "treverter2d.hpp"
#include "healgrid.hpp"
#include "buildgrid.hpp"
#include "modgrid.hpp"
#include "finder2d.hpp"
#include "Gmsh.h"
#include "GModel.h"
#include "MVertex.h"
#include "MPoint.h"
#include "MLine.h"
#include "nan_handler.h"
using namespace HM2D;

HMCallback::FunctionWithCallback<Mesher::TUnstructuredTriangle> Mesher::UnstructuredTriangle;
HMCallback::FunctionWithCallback<Mesher::TUnstructuredTriangleRecomb> Mesher::UnstructuredTriangleRecomb;

namespace{
Contour::Tree no_crosses_with_priority(const Contour::Tree& source){
	Contour::Tree ret = Contour::Tree::DeepCopy(source, 2);

	//points which should be kept in result tree
	std::vector<Point*> priority;
	for (auto& p: AllVertices(source.alledges())) priority.push_back(p.get());
	std::sort(priority.begin(), priority.end());

	auto get_priority = [&priority](const std::tuple<bool, shared_ptr<Vertex>>& res)->int{
		if (std::get<0>(res) == false) return 0;

		Point* p = std::get<1>(res).get();
		if (std::binary_search(priority.begin(), priority.end(), p))
			return 1;

		return 0;
	};

	//add all crosses points
	auto treat_conts = [&get_priority](HM2D::EdgeData& c1, HM2D::EdgeData& c2){
		auto crosses = HM2D::Contour::Finder::CrossAll(c1, c2);
		VertexData from1, from2, to;
		for (auto& c: crosses){
			auto res1 = HM2D::Contour::GuaranteePoint(c1, std::get<1>(c));
			auto res2 = HM2D::Contour::GuaranteePoint(c2, std::get<1>(c));
			if (std::get<1>(res1) == std::get<1>(res2)) continue;
			//substitute res1, res2 with node of highest priority
			int pri1 = get_priority(res1);
			int pri2 = get_priority(res2);
			shared_ptr<Vertex> hp = std::get<1>(res1);
			if (pri2 > pri1) hp = std::get<1>(res2);
			from1.push_back(std::get<1>(res1));
			from2.push_back(std::get<1>(res2));
			to.push_back(hp);
		}
		if (from1.size() == 0) return;
		//do substitutions
		aa::constant_ids_pvec(HM2D::AllVertices(c1), -1);
		aa::constant_ids_pvec(HM2D::AllVertices(c2), -1);
		aa::enumerate_ids_pvec(from1);
		aa::enumerate_ids_pvec(from2);
		aa::constant_ids_pvec(to, -1);
		for (auto e: c1){
			if (e->first()->id != -1) e->vertices[0] = to[e->first()->id];
			if (e->last()->id != -1) e->vertices[1] = to[e->last()->id];
		}
		for (auto e: c2){
			if (e->first()->id != -1) e->vertices[0] = to[e->first()->id];
			if (e->last()->id != -1) e->vertices[1] = to[e->last()->id];
		}
	};
	for (int i=0; i<ret.nodes.size(); ++i)
	for (int j=0; j<i; ++j){
		treat_conts(ret.nodes[i]->contour, ret.nodes[j]->contour);
	}

	return ret;
}
}

Contour::Tree Mesher::PrepareSource(const Contour::Tree& source, double defsize){
	//paste crosses
	Contour::Tree ret = no_crosses_with_priority(source);

	//sizes calculation
	std::map<Point*, double> sizes;
	auto ae = ret.alledges();
	auto lens = ELengths(ae);
	for (auto ve: Connectivity::VertexEdge(ae)){
		double lenmin=lens[ve.eind[0]], lenmax=lens[ve.eind[0]];
		for (int i=1; i<ve.size(); ++i){
			double len = lens[ve.eind[i]];
			if (len < lenmin) lenmin = len;
			if (len > lenmax) lenmax = len;
		}
		sizes[ve.v.get()] = std::min(lenmax, 2.5*lenmin);
	}

	//defsize
	if (defsize>0) for (auto& it: sizes){
		if (it.second > defsize) it.second = defsize;
	}
	//do partitions
	for (auto& n: ret.nodes){
		auto w = HM2D::Contour::EWeights(n->contour);
		auto op = HM2D::Contour::OrderedPoints(n->contour);
		std::map<double, double> basis;
		for (int i=0; i<w.size(); ++i){
			basis[w[i]] = sizes[op[i].get()];
		}
		//partition
		EdgeData newcont = Contour::Algos::WeightedPartition(basis, n->contour,
				Contour::Algos::PartitionTp::KEEP_ALL);
		std::swap(n->contour, newcont);
	}

	return ret;
}
Contour::Tree Mesher::PrepareSource(const EdgeData& source, double defsize){
	Contour::Tree tree = Contour::Tree::Assemble(source);
	return PrepareSource(tree, defsize);
}

EdgeData Mesher::RepartSourceById(const EdgeData& source){
	EdgeData ret;
	auto it = source.begin();
	VertexData keep_pts;
	for (auto v: AllVertices(source)) if (v->id == 1){
		keep_pts.push_back(v);
	}
	while (it != source.end()){
		if ((*it)->id == 1){
			ret.push_back(std::make_shared<Edge>(**it));
			++it;
		} else {
			auto itprev = it;
			if (it == source.begin()) itprev = std::prev(source.end());
			else --itprev;
			EdgeData group;
			do{
				group.push_back(*it++);
			} while (it != source.end() && (*it)->id != 1);
			auto itnext = it;
			if (itnext == source.end()) itnext = source.begin();
			std::map<double, double> weights;
			weights.emplace(0, (*itprev)->length());
			weights.emplace(1, (*itnext)->length());
			EdgeData repart = Contour::Algos::WeightedPartition(weights, group,
					keep_pts);
			ret.insert(ret.end(), repart.begin(), repart.end());
		}
	}

	return ret;
}

namespace{

GridData GridFromModel(GModel& m){
	m.indexMeshVertices(true);
	VertexData vvert(m.getNumMeshVertices());

	vector<vector<int>> cellvert;
	for (auto it = m.firstFace(); it!=m.lastFace(); ++it){
		for (int en=0; en<(*it)->getNumMeshElements(); ++en){
			auto e = (*it)->getMeshElement(en);
			std::vector<MVertex*> elvert;
			e->getVertices(elvert);
			cellvert.emplace_back(elvert.size());
			for (int i=0; i<e->getNumVertices(); ++i){
				auto v = elvert[i];
				int index = v->getIndex()-1;
				assert(index<vvert.size());
				cellvert.back()[i] = index;
				if (!vvert[index]) vvert[index].reset(new Vertex(v->x(), v->y()));
			}
		}
	}

	return Grid::Constructor::FromTab(std::move(vvert), cellvert);
}

vector<Contour::Tree> build_cropped(const Contour::Tree& source){
	//using shallow copies
	vector<Contour::Tree> ret;
	//constraints preproc
	vector<BoundingBox> constrbb;
	for (auto& n: source.detached_contours()){
		constrbb.push_back(HM2D::BBox(n->contour));
	}
	for (auto& n: source.nodes) if (n->level % 2 == 0){
		//root
		ret.emplace_back();
		auto& t = ret.back();
		t.nodes.push_back(n);
		//children
		for (auto& ch: n->children){
			t.nodes.push_back(ch.lock());
		}
		//constraints for this subtree: quick bounding box check only
		if (constrbb.size()==0) continue;
		
		vector<BoundingBox> nbb;
		for (int i=0; i<t.nodes.size(); ++i){
			nbb.push_back(BBox(n->contour));
		}

		int ic = 0;
		for (auto& n: source.detached_contours()){
			auto& cbb = constrbb[ic++];
			int rel = cbb.relation(nbb[0]);
			if (rel == 3) continue;

			for (int i=1; i<t.nodes.size(); ++i){
				rel = cbb.relation(nbb[i]);
				if (rel == 2) break;
			}
			if (rel == 2) continue;

			t.nodes.push_back(n);
		}
	}

	return ret;
}

vector<GEdge*> fill_model_with_1d(
		GModel& m,
		const EdgeData& edata,
		vector<GEdge*>& eheap,
		vector<GVertex*>& vheap,
		vector<MVertex*>& mvheap,
		double dist){
	auto place_vertex = [&](Vertex& v)->GVertex*{
		if (vheap[v.id] == 0){
			vheap[v.id] = m.addVertex(v.x, v.y, 0, dist);
			mvheap[v.id] = new MVertex(v.x, v.y, 0, vheap[v.id]);
			vheap[v.id]->mesh_vertices.push_back(mvheap[v.id]);
			vheap[v.id]->points.push_back(new MPoint(mvheap[v.id]));
		}
		auto ms = vheap[v.id]->prescribedMeshSizeAtVertex();
		if (ms>dist) vheap[v.id]->setPrescribedMeshSizeAtVertex(dist);
		return vheap[v.id];
	};
	auto place_edge = [&](Edge& ecur, Edge& eprev)->GEdge*{
		if (eheap[ecur.id] != 0) return eheap[ecur.id];
		auto v1 = ecur.first().get();
		auto v2 = ecur.last().get();
		if (v1 != eprev.first().get() && v1 != eprev.last().get()){
			std::swap(v1, v2);
		}
		auto gv1 = place_vertex(*v1);
		auto gv2 = place_vertex(*v2);
		eheap[ecur.id] = m.addLine(gv1, gv2);
		eheap[ecur.id]->addLine(new MLine(mvheap[v1->id], mvheap[v2->id]));
		return eheap[ecur.id];
	};
	vector<GEdge*> ret;
	Edge* prev = edata.back().get();
	for (auto& e: edata){
		ret.push_back(place_edge(*e, *prev));
		prev = e.get();
	}
	return ret;
}

GFace* assemble_face(GModel& m, const Contour::Tree& tree, vector<vector<GEdge*>>& g_edges,
		const std::map<Point, double>& embedded){
	vector<vector<GEdge*>> fedges(1);
	//root
	for (int i=0; i<tree.nodes.size(); ++i) if (tree.nodes[i]->level % 2 == 0){
		std::swap(fedges[0], g_edges[i]);
		break;
	}
	//children
	for (int i=0; i<tree.nodes.size(); ++i) if (tree.nodes[i]->level % 2 == 1){
		fedges.emplace_back();
		std::swap(fedges.back(), g_edges[i]);
	}
	GFace* fc = m.addPlanarFace(fedges);

	//constraints
	for (int i=0; i<tree.nodes.size(); ++i) if (tree.nodes[i]->isdetached()){
		for (auto e: g_edges[i]) fc->addEmbeddedEdge(e);
	}
	//embedded points
	for (auto& p: embedded){
		if (tree.whereis(p.first) != INSIDE) continue;
		auto added = m.addVertex(p.first.x, p.first.y, 0, p.second);
		fc->addEmbeddedVertex(added);
	}

	return fc;
}

void fill_model_with_2d(GModel& m, const GFace* fc){
	//Mesh2D
	//m.writeGEO("gmsh_geo.geo");
	//m.writeMSH("gmsh_geo.msh");
	m.mesh(2);
	//m.writeVTK("gmsh_geo.vtk");
	//m.writeMSH("gmsh_geo.msh");
}

void fill_model_with_2d_recomb(GModel& m, GFace* fc){
	//usage of delaunay for quads gives worse results for non-regular areas
	//hence using auto algorithm
	GmshSetOption("Mesh", "Algorithm", 2.0);
	GmshSetOption("Mesh", "RecombinationAlgorithm", 1.0);
	fc->meshAttributes.recombine = 1.0;
	m.mesh(2);

	// if all nodes of quad grid are still trianlge than most likely builder
	// has failed. So we use another algorithm
	bool has4=false;
	for (int en=0; en<fc->getNumMeshElements(); ++en){
		auto e = fc->getMeshElement(en);
		if (e->getNumVertices() == 4) { has4=true; break; }
	}
	if (!has4){
		GmshSetOption("Mesh", "Algorithm", 2.0);
		GmshSetOption("Mesh", "RecombinationAlgorithm", 0.0);
		fc->meshAttributes.recombine = 1.0;
		m.mesh(2);
	}
}

GridData gmsh_fill(const Contour::Tree& tree, const std::map<Point, double>& embedded,
		int algo, HMCallback::Caller2& cb){
	//2 - auto, 5 - delaunay, 6 - frontal, 8 - delaunay for quads
	GmshSetOption("Mesh", "Algorithm", 2.0);
	GmshSetOption("Mesh", "Optimize", 1.0);
	GmshSetOption("General", "Verbosity", 0.0);
	GModel m;
	m.setFactory("Gmsh");

	auto ae = tree.alledges();
	auto av = AllVertices(ae);
	aa::enumerate_ids_pvec(ae);
	aa::enumerate_ids_pvec(av);

	//size
	double h = 2*HM2D::BBox(tree.alledges()).lendiag();

	cb.step_after(10, "Fill 1D mesh");
	vector<GEdge*> g_edges_heap(ae.size(), 0);
	vector<GVertex*> g_vertex_heap(av.size(), 0);
	vector<MVertex*> m_vertex_heap(av.size(), 0);
	vector<vector<GEdge*>> g_edges;
	for (auto n: tree.nodes){
		g_edges.push_back(fill_model_with_1d(m, n->contour, g_edges_heap, g_vertex_heap, m_vertex_heap, h));
	}
	cb.step_after(10, "Face assembling");
	GFace* gf = assemble_face(m, tree, g_edges, embedded);

	cb.step_after(50, "Meshing");
	//!! gmsh 2.11 doesn't work correctly without this line
	//   if chararcteristic mesh size is not 1.0
	auto bb = m.bounds();
	GmshSetBoundingBox(bb.min()[0], bb.max()[0], bb.min()[1], bb.max()[1], 0, 0);
	NanSignalHandler::StopCheck();
	if (algo == 0) fill_model_with_2d(m, gf);
	else fill_model_with_2d_recomb(m, gf);
	NanSignalHandler::StartCheck();

	cb.step_after(30, "Assemble mesh");
	return GridFromModel(m);
}

//uses 100 units of callback
GridData gmsh_builder(const Contour::Tree& source, const std::map<Point, double>& embedded, int algo,
		shared_ptr<HMCallback::Caller2> callback){
	if (source.roots().size() == 0) return GridData(); 

	callback->step_after(10, "Prepare contours");
	HM2D::Contour::R::RevertTree rt(source);

	vector<Contour::Tree> trees = build_cropped(source);
	GridData ret;

	for (int i=0; i<trees.size(); ++i){
		auto cb = callback->subrange(85./trees.size(), 100.);
		if (i == 0) ret = gmsh_fill(trees[i], embedded, algo, *cb);
		else {
			GridData sg = gmsh_fill(trees[i], embedded, algo, *cb);
			std::copy(sg.vvert.begin(), sg.vvert.end(), std::back_inserter(ret.vvert));
			std::copy(sg.vedges.begin(), sg.vedges.end(), std::back_inserter(ret.vedges));
			std::copy(sg.vcells.begin(), sg.vcells.end(), std::back_inserter(ret.vcells));
		}
	}

	callback->step_after(5, "Finalizing");
	//adjust rotation
	for (int i=0; i<ret.vcells.size(); ++i){
		if (Contour::Area(ret.vcells[i]->edges) < 0){
			Contour::Reverse(ret.vcells[i]->edges);
		}
	}
	return ret;
}

}

GridData Mesher::TUnstructuredTriangle::_run(const Contour::Tree& source, const std::map<Point, double>& embedded){
	return gmsh_builder(source, embedded, 0, callback);
}

GridData Mesher::TUnstructuredTriangle::_run(const Contour::Tree& source){
	return _run(source, std::map<Point, double>());
}

namespace{

void recomb_heal(GridData& grid){
	auto process = [&](int ic1, int ic2, int n1, int n2){
		Cell* c1 = grid.vcells[ic1].get();
		Cell* c2 = grid.vcells[ic2].get();
		Contour::R::ReallyDirect::Permanent(c1->edges);
		Contour::R::ReallyRevert::Permanent(c2->edges);
		int loc1 = -1, loc2 = -1;
		for (int j=0; j<c1->edges.size(); ++j){
			if (c1->edges[j]->vertices[0]->id == n1){
				loc1 = j; break;
			}
		}
		for (int j=0; j<c2->edges.size(); ++j){
			if (c2->edges[j]->vertices[0]->id == n1){
				loc2 = j; break;
			}
		}
		assert(loc1 >= 0 && loc2 >= 0);
		assert(c1->edges[loc1]->vertices[1]->id == n2);
		std::rotate(c2->edges.begin(), c2->edges.begin()+loc2, c2->edges.end());
		c1->edges.insert(c1->edges.begin()+loc1+1, c2->edges.begin()+1, c2->edges.end()-1);
		grid.vcells[ic2]=nullptr;
		Grid::Algos::RestoreFromCells(grid);
		return !std::get<0>(Contour::Finder::SelfCross(c1->edges));
	};

	for (int tries=0; tries<100; ++tries){
		vector<std::map<int, int>> used_edges(grid.vvert.size());
		aa::enumerate_ids_pvec(grid.vvert);
		for (int ic=0; ic<grid.vcells.size(); ++ic){
			Cell* c=grid.vcells[ic].get();
			auto op = Contour::OrderedPoints(c->edges);
			for (int i=0; i<op.size()-1; ++i){
				int p1 = op[i]->id;
				int p2 = op[i+1]->id;
				auto er = used_edges[p1].emplace(p2, ic);
				if (er.second==false){
					if (process(ic, er.first->second, p1, p2)) goto NEXT_TRY;
					else goto ERR_OUT;
				}
			}
		}
		//no changes in current try => grid is fine.
		if (tries > 0) Grid::Algos::CutCellDims(grid, 4);
		return;
NEXT_TRY:
		continue; 
	}
ERR_OUT:
	throw std::runtime_error("cannot restore correct grid from gmsh output");

}

void guarantee_edges(GridData& grid, const EdgeData& edges){
	if (edges.size() == 0) return;
	//find all edges vertices in the grid
	auto av = AllVertices(edges);
	VertexData vfnd;
	Finder::VertexMatch vfinder(grid.vvert);
	for (auto& v: av) vfnd.push_back(vfinder.find(*v));
	if (std::find(vfnd.begin(), vfnd.end(), nullptr) != vfnd.end()){
		throw std::runtime_error("failed to satisfy constraint condition");
	}

	//find all constraint edges in the grid
	//if fail -> write node pairs for splitting
	vector<std::pair<Vertex*, Vertex*>> force_edges;
	Finder::EdgeFinder efinder(grid.vedges);
	aa::enumerate_ids_pvec(av);
	for (auto e: edges){
		Vertex* gv1 = vfnd[e->vertices[0]->id].get();
		Vertex* gv2 = vfnd[e->vertices[1]->id].get();
		auto ef = efinder.find(gv1, gv2);
		if (!std::get<0>(ef)) force_edges.emplace_back(gv1, gv2);
	}
	if (force_edges.size() == 0) return;

	//get cell index for each split pair
	//and local indicies
	vector<int> icells, lnode1, lnode2;
	aa::constant_ids_pvec(grid.vvert, -1);
	for (int i=0; i<force_edges.size(); ++i){
		force_edges[i].first->id = i;
		force_edges[i].second->id = i;
	}
	for (int ic=0; ic<grid.vcells.size(); ++ic) if (grid.vcells[ic]->edges.size()>3){
		auto op = Contour::OrderedPoints(grid.vcells[ic]->edges);
		for (int i=0; i<op.size()-3; ++i){
			if (op[i]->id >=0){
				for (int j=i+2; j<op.size()-1; ++j){
					if (op[j]->id == op[i]->id){
						icells.push_back(ic);
						lnode1.push_back(i);
						lnode2.push_back(j);
						break;
					}
				}
				break;
			}
		}
		if (icells.size() == force_edges.size()) break;
	}

	//split unique icells
	std::set<int> unique_icells;
	for (int i=0; i<icells.size(); ++i){
		if (unique_icells.insert(icells[i]).second){
			Grid::Algos::SplitCell(grid, icells[i], lnode1[i], lnode2[i]);
		}
	}
}

};

GridData Mesher::TUnstructuredTriangleRecomb::_run(const Contour::Tree& source, const std::map<Point, double>& embedded){
	GridData ret = gmsh_builder(source, embedded, 1, callback);

	callback->step_after(10, "Recombination check");
	//as a result of recombination some narrow reversed boundary triangles may occur.
	//here we try to fix it by merging with adjacent inner cells.
	recomb_heal(ret);
	//now we should check constraints lay on grid edges
	//since this feature can be not satisfied after recombination.
	EdgeData ec;
	for (auto n: source.nodes) if (n->isdetached()){
		ec.insert(ec.end(), n->contour.begin(), n->contour.end());
	}
	guarantee_edges(ret, ec);
	
	return ret;
}

GridData Mesher::TUnstructuredTriangleRecomb::_run(const Contour::Tree& source){
	return _run(source, std::map<Point, double>());
}
