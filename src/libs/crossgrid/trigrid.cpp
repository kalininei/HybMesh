#include "trigrid.hpp"
#include "modcont.hpp"
#include "partcont.hpp"
#include "contabs2d.hpp"
#include "treverter2d.hpp"
#include "healgrid.hpp"
#include "buildgrid.hpp"
#include "modgrid.hpp"
#include "finder2d.hpp"
#include "buildcont.hpp"
#include "assemble2d.hpp"
#include "Gmsh.h"
#include "GModel.h"
#include "MVertex.h"
#include "MPoint.h"
#include "MLine.h"
#include "nan_handler.h"
#include "GmshMessage.h"
using namespace HM2D;

HMCallback::FunctionWithCallback<Mesher::TUnstructuredTriangle> Mesher::UnstructuredTriangle;
HMCallback::FunctionWithCallback<Mesher::TUnstructuredTriangleRecomb> Mesher::UnstructuredTriangleRecomb;

namespace{

class GmshCallback: public GmshMessage{
	void operator()(std::string level, std::string message) override{
		std::cout<<level<<": "<<message<<std::endl;
	}
};
GmshCallback gmshcb;

}

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
			auto res1 = HM2D::Contour::Algos::GuaranteePoint(c1, std::get<1>(c));
			auto res2 = HM2D::Contour::Algos::GuaranteePoint(c2, std::get<1>(c));
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
	auto lens = Contour::ELengths(ae);
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

namespace{

std::map<double, double> build_map01(double sz0, double sz1){
	assert(sz0>0 && sz1>0);
	std::map<double, double> ret;
	ret[0] = sz0; ret[1] = sz1;
	return ret;
}

std::map<double, double> fill_map_keys(double contlen, double szend,
		const vector<std::pair<Point, double>>& src, int nmax){
	double minsz;
	if (szend>=0) minsz=szend;
	else minsz = src[0].second;
	for (int i=0; i<src.size(); ++i){
		if (src[i].second < minsz) minsz = src[i].second;
	}
	int n = 0.5 * contlen / minsz;
	if (n < 5) n = 5;
	if (n > nmax) n = nmax;
	std::map<double, double> ret;
	for (int i=0; i<n+1; ++i){
		ret.emplace(double(i)/n,  0);
	}
	return ret;
}

bool not_intersect(Point p1, Point p2, const EdgeData& cont){
	p1 = Point::Weigh(p1, p2, 0.001);
	p2 = Point::Weigh(p1, p2, 0.999);
	auto c2 = Contour::Constructor::FromPoints({p1, p2});
	return !std::get<0>(Contour::Finder::Cross(cont, c2));
}

bool is_valueble(Point p1, const EdgeData& cont){
	auto c = Contour::CoordAt(cont, p1);
	return ISIN_NN(std::get<1>(c), 0, 1);
}

void find_closest(const EdgeData& cont, const Point& p, const vector<std::pair<Point, double>>& src,
		double& dist, double& sz){
	dist = 1e200; sz = 1e200; Point pbest;
	auto ce = Finder::ClosestEdge(cont, p);
	Point *p1 = cont[std::get<0>(ce)]->pfirst();
	Point *p2 = cont[std::get<0>(ce)]->plast();
	if (!Contour::CorrectlyDirectedEdge(cont, std::get<0>(ce))){
		std::swap(p1, p2);
	}
	for (auto& s: src){
		if (LinePointWhereIs(s.first, *p1, *p2) == 0){
			double d = Point::meas(p, s.first);
			if (d < dist &&
			    not_intersect(p, s.first, cont)){
				pbest = s.first;
				dist = d;
				sz = s.second;
			}
		}
	}

	if (dist != 1e200 && is_valueble(pbest, cont)){
		dist = sqrt(dist);
	} else { dist = -1; sz = -1; }
}

double amean(double l0, double sz0, double l1, double sz1){
	return (l0*sz1 + l1*sz0)/(l0+l1);
}
double amean(double l0, double sz0, double l1, double sz1, double l2, double sz2){
	return ((l1+l2)*sz0 + (l0+l2)*sz1 + (l0+l1)*sz2)/(l0+l1+l2)/2.;
}
double calculate_len_at(const EdgeData& cont, const Point& p, double l0, double sz0, double l1, double sz1,
		const vector<std::pair<Point, double>>& src){
	if (ISZERO(l0)) return sz0;
	if (ISZERO(l1)) return sz1;
	double d, s;
	find_closest(cont, p, src, d, s);
	if (d<0) return amean(l0, sz0, l1, sz1);
	if (ISZERO(d)) return d;
	//if one point is far away
	if (d>l0+l1) return amean(l0, sz0, l1, sz1);
	//if (d<l0 && l0<l1) return amean(l0, sz0, d, s);
	//if (d<l1 && l1<l0) return amean(l1, sz1, d, s);
	//calc weights
	return amean(l0, sz0, l1, sz1, d, s);
}

double calculate_len_at(const EdgeData& cont, const Point& p, const vector<std::pair<Point, double>>& src){
	double d, s;
	find_closest(cont, p, src, d, s);
	return s;
}

std::map<double, double> build_weights(const EdgeData& cont, double sz0, double sz1,
		const vector<std::pair<Point, double>>& src){
	if (src.size() == 0) return build_map01(sz0, sz1);
	double len = HM2D::Contour::Length(cont);
	auto ret = fill_map_keys(len, std::min(sz0, sz1), src, 100);

	vector<double> w;
	for (auto& r: ret) w.push_back(r.first);
	vector<Point> points = Contour::WeightPoints(cont, w);
	if (sz0 >=0){
		ret.begin()->second = sz0;
		ret.rbegin()->second = sz1;
		auto it = std::next(ret.begin());
		for (int i=1; i<ret.size()-1; ++i, ++it){
			it->second = calculate_len_at(cont, points[i],
				it->first*len, sz0, (1-it->first)*len, sz1,
				src);
		}
	} else {
		auto it = ret.begin();
		for (int i=0; i<points.size(); ++i){
			it->second = calculate_len_at(cont, points[i], src);
			if ((*it++).second < 0) ret.erase(std::prev(it));
		}
	}
	return ret;
};

EdgeData repart_cont(const EdgeData& cont, double sz0, double sz1, 
		const VertexData& keep_pts,
		const vector<std::pair<Point, double>>& src){
	std::map<double, double> weights = build_weights(cont, sz0, sz1, src);
	return Contour::Algos::WeightedPartition(weights, cont, keep_pts);
}

vector<std::pair<Point, double>> sort_out_sources(const EdgeData& source,
		const vector<std::pair<Point, double>>& input){
	aa::RestoreIds<EdgeData> rr(source);

	vector<Point> pp(input.size());
	for (int i=0; i<input.size(); ++i) pp[i] = input[i].first;
	vector<int> srt = Contour::Finder::SortOutPoints(source, pp);
	vector<std::pair<Point, double>> ret;
	for (int i=0; i<srt.size(); ++i) if (srt[i] != OUTSIDE){
		ret.push_back(input[i]);
	}
	return ret;
}

VertexData assemble_keep_pts(const EdgeData& source){
	VertexData keep_pts;
	for (auto e: source) if (e->id != 1){
		if (e->pfirst()->id == 1) {
			keep_pts.push_back(e->first());
			e->pfirst()->id = 0;
		}
		if (e->plast()->id == 1) {
			keep_pts.push_back(e->last());
			e->plast()->id = 0;
		}
	}
	return keep_pts;
}

EdgeData start_from_id1(const EdgeData& source){
	EdgeData src2 = source;
	for (auto it = src2.begin(); it != src2.end(); ++it){
		if ((*it)->id == 1){
			std::rotate(src2.begin(), it, src2.end());
			break;
		}
	}
	return src2;
}
};

EdgeData Mesher::RepartSourceById(const EdgeData& source,
		const vector<std::pair<Point, double>>& size_src){
	assert(Contour::IsClosed(source));
	Contour::R::Clockwise cc(source, false);
	EdgeData ret;
	//if nothing to part return deepcopy
	if (std::all_of(source.begin(), source.end(), [](const shared_ptr<Edge>& e){
			return e->id == 1; })){
		DeepCopy(source, ret, 0);
		return ret;
	}
	//keep pts
	VertexData keep_pts = assemble_keep_pts(source);

	//sort out only those size_src which lie inside source
	vector<std::pair<Point, double>>
	size_src2 = sort_out_sources(source, size_src);

	//find any edge with id = 1 to start assembling from it
	EdgeData src2 = start_from_id1(source);

	//start partition
	auto it = src2.begin();
	while (it != src2.end()){
		if ((*it)->id == 1){
			ret.push_back(std::make_shared<Edge>(**it));
			++it;
		} else {
			auto itprev = it;
			EdgeData group;
			do{
				group.push_back(*it++);
			} while (it != src2.end() && (*it)->id != 1);
			//start/end sizes
			double sz0=-1, sz1=-1;
			if (itprev != src2.begin()) sz0 = (*(itprev-1))->length();
			if (it != src2.end()) sz1 = (*it)->length();
			else sz1 = (*src2.begin())->length();
			//repartition
			EdgeData repart = repart_cont(group, sz0, sz1, keep_pts, size_src2);
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

namespace{
struct ShpVCmp{
	bool operator()(const shared_ptr<HM2D::Vertex>& v1, const shared_ptr<HM2D::Vertex>& v2){
		return *v1 < *v2;
	}
};
};

vector<Contour::Tree> build_cropped(const Contour::Tree& source){
	//using shallow copies
	vector<Contour::Tree> ret;
	//constraints preproc
	vector<BoundingBox> constrbb;
	for (auto& n: source.detached_contours()){
		constrbb.push_back(HM2D::BBox(n->contour));
	}
	for (auto& n: source.nodes) if (n->isouter()){
		//root
		Contour::Tree t;
		t.nodes.push_back(n);
		//children
		for (auto& ch: n->children){
			t.nodes.push_back(ch.lock());
		}
		//constraints for this subtree: quick bounding box check only
		if (constrbb.size()!=0){
			vector<BoundingBox> nbb;
			for (int i=0; i<t.nodes.size(); ++i){
				nbb.push_back(BBox(t.nodes[i]->contour));
			}
			int ic = 0;
			for (auto& n: source.detached_contours()){
				auto& cbb = constrbb[ic++];
				if (cbb.relation(nbb[0]) == 3) continue;
				t.nodes.push_back(n);
			}
		}

		//merge equal points
		std::set<shared_ptr<Vertex>, ShpVCmp> vset;
		t = Contour::Tree::DeepCopy(t);
		auto merge_ep = [&vset](shared_ptr<HM2D::Vertex>& v){
			auto er = vset.insert(v);
			if (!er.second){ v = *er.first; }
		};
		for (auto& e: t.alledges()){
			merge_ep(e->vertices[0]);
			merge_ep(e->vertices[1]);
		}

		ret.push_back(std::move(t));
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
		const CoordinateMap2D<double>& embedded){
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
	for (auto it=embedded.begin(); it!=embedded.end(); ++it){
		if (tree.whereis(it.point()) != INSIDE) continue;
		auto added = m.addVertex(it.x(), it.y(), 0, it.data());
		auto mv = new MVertex(it.x(), it.y(), 0, added);
		added->mesh_vertices.push_back(mv);
		added->points.push_back(new MPoint(mv));
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

	// if all nodes of quad grid are still trianlge then most likely builder
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

GridData gmsh_fill(const Contour::Tree& tree, const CoordinateMap2D<double>& embedded,
		int algo, HMCallback::Caller2& cb){
	GModel m;
	m.setFactory("Gmsh");
	//2 - auto, 5 - delaunay, 6 - frontal, 8 - delaunay for quads
	GmshSetOption("Mesh", "Algorithm", 2.0);
	GmshSetOption("Mesh", "Optimize", 1.0);
	GmshSetOption("General", "Verbosity", 0.0);
	//GmshSetOption("General", "Verbosity", 100.0);
	//GmshSetMessageHandler(&gmshcb);

	auto ae = tree.alledges();
	auto av = AllVertices(ae);
	//enumeration will be used in fill_model_with_1d procedure
	aa::enumerate_ids_pvec(ae);
	aa::enumerate_ids_pvec(av);

	//to prevent creating new boundary nodes
	//default size passed to gmsh is bigger than any edge length.
	//not sure if i need it since i explicitly define 1d mesh.
	double h = 2*HM2D::BBox(av).lendiag();

	cb.step_after(10, "Fill 1D mesh");
	vector<GEdge*> g_edges_heap(ae.size(), 0);
	vector<GVertex*> g_vertex_heap(av.size(), 0);
	vector<MVertex*> m_vertex_heap(av.size(), 0);
	vector<vector<GEdge*>> g_edges;


	for (auto n: tree.nodes){
		g_edges.push_back(fill_model_with_1d(
			m, n->contour, g_edges_heap, g_vertex_heap, m_vertex_heap, h));
	}
	cb.step_after(10, "Face assembling");
	GFace* gf = assemble_face(m, tree, g_edges, embedded);
	

	cb.step_after(50, "Meshing");
	//!! gmsh 2.11 doesn't work correctly without this line
	//   if chararcteristic mesh size is not 1.0
	auto bb = m.bounds();
	GmshSetBoundingBox(bb.min()[0], bb.max()[0], bb.min()[1], bb.max()[1], 0, 0);
	//gmsh has some zero division operations hence we need to stop checking
	NanSignalHandler::StopCheck();
	if (algo == 0) fill_model_with_2d(m, gf);
	else fill_model_with_2d_recomb(m, gf);
	//turn nan check on
	NanSignalHandler::StartCheck();

	cb.step_after(25, "Assemble mesh");
	GridData ret = GridFromModel(m);

	cb.step_after(5, "Assign boundary types");
	EdgeData gbnd = HM2D::ECol::Assembler::GridBoundary(ret);
	EdgeData cbnd;
	for (auto& n: tree.bound_contours())
		cbnd.insert(cbnd.end(), n->contour.begin(), n->contour.end());
	HM2D::ECol::Algos::AssignBTypes(cbnd, gbnd);

	return ret;
}

//uses 100 units of callback
GridData gmsh_builder(const Contour::Tree& source, const CoordinateMap2D<double>& embedded, int algo,
		shared_ptr<HMCallback::Caller2> callback){
	if (source.roots().size() == 0) return GridData(); 

	callback->step_after(10, "Prepare contours");
	HM2D::Contour::R::RevertTree rt(source);

	vector<Contour::Tree> trees = build_cropped(source);

	vector<GridData> gg;
	for (int i=0; i<trees.size(); ++i){
		auto cb = callback->subrange(80./trees.size(), 100.);
		gg.push_back(gmsh_fill(trees[i], embedded, algo, *cb));
	}

	callback->step_after(10, "Finalizing");
	//adjust rotation
	for (auto ret: gg)
	for (int i=0; i<ret.vcells.size(); ++i){
		if (Contour::Area(ret.vcells[i]->edges) < 0){
			Contour::Algos::Reverse(ret.vcells[i]->edges);
		}
	}
	GridData ret = std::move(gg[0]);
	//make a summation with merging
	for (int i=1; i<gg.size(); ++i){
		Grid::Algos::MergeBoundaries(gg[i], ret);
	}
	//merge equal boundary points
	return ret;
}

}

GridData Mesher::TUnstructuredTriangle::_run(const Contour::Tree& source, const CoordinateMap2D<double>& embedded){
	return gmsh_builder(source, embedded, 0, callback);
}

GridData Mesher::TUnstructuredTriangle::_run(const Contour::Tree& source){
	return _run(source, CoordinateMap2D<double>());
}

namespace{
//This is a heavy workaround for gmsh recombination algorithm failures.
void recomb_heal(GridData& grid){
	auto get_adjacent_cells = [&](Edge* e, int& c1, int& c2){
		aa::enumerate_ids_pvec(grid.vcells);
		c1 = -1; c2 = -1;
		if (e->has_left_cell()) c1=e->left.lock()->id;
		if (e->has_right_cell()) c2=e->right.lock()->id;
		//we cannot rely of left/right only because
		//self crossing cells with unspecified edge direction can be passed here.
		if (c2 == c1) c2 = -1;
		if ((c1>-1 && c2>-1) || (c1==-1 && c2==-1)) return;
		c1 = std::max(c1, c2); c2 = -1;
		for (int i=0; i<grid.vcells.size(); ++i){
			if (HM2D::Finder::Contains(grid.vcells[i]->edges, e) && i!=c1){
				c2 = i;
				return;
			}
		}
	};
	std::function<bool(int, int)> process = [&](int ic1, int ic2)->bool{
		//This procedire connects ic1 and ic2 cells by removing common
		//edges and reassembling other edges.
		Cell* c1 = grid.vcells[ic1].get();
		Cell* c2 = grid.vcells[ic2].get();
		aa::constant_ids_pvec(c1->edges, 0);
		aa::constant_ids_pvec(c2->edges, 0);
		for (auto e: c1->edges) ++e->id;
		for (auto e: c2->edges) ++e->id;
		EdgeData ee;
		for (auto e: c1->edges) if (e->id == 1) ee.push_back(e);
		for (auto e: c2->edges) if (e->id == 1) ee.push_back(e);
		ee = HM2D::Contour::Assembler::Contour1(ee);
		//if ee is not closed cell return error status.
		if (!HM2D::Contour::IsClosed(ee)) return false;
		Contour::R::Clockwise::Permanent(ee, false);
		c1->edges = ee;
		for (auto e: c1->edges){
			e->left = grid.vcells[ic1];
			if (e->right.lock().get() == c2) e->right.reset();
		}
		grid.vcells[ic2]=nullptr;
		Grid::Algos::RestoreFromCells(grid);
		auto cc = Contour::Finder::SelfCross(ee);
		// If resuling cell has no self crosses it's ok to quit the function.
		// See others/botscript1.py
		if (!std::get<0>(cc)) return true;
		// if there is a cross we have to analyse each of crossed edges
		// and delete it if it has two adjacent cells.
		// See others/unite_grids1.py, others/rrj2 (UNITE1 and UNITE4).
		get_adjacent_cells(ee[std::get<4>(cc)].get(), ic1, ic2);
		if (ic2 > -1) return process(ic1, ic2);
		get_adjacent_cells(ee[std::get<5>(cc)].get(), ic1, ic2);
		if (ic2 > -1) return process(ic1, ic2);
		return false;
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
					if (process(ic, er.first->second)) goto NEXT_TRY;
					else goto ERR_OUT;
				}
			}
		}
		//no changes in current try => grid is fine.
		if (tries > 0){
			Grid::Algos::CutCellDims(grid, 4);
			Grid::Algos::NoConcaveCells(grid, 175);
		}
		return;
NEXT_TRY:
		continue; 
	}
ERR_OUT:
	throw std::runtime_error("cannot restore correct grid from gmsh recombination procedure output");
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

	//remove all cell edges which cross force_edges
	aa::constant_ids_pvec(grid.vvert, -1);
	for (int i=0; i<force_edges.size(); ++i){
		force_edges[i].first->id = i;
		force_edges[i].second->id = i;
	}
	vector<int> bad_edges;
	auto bb = HM2D::BBox(grid.vvert);
	BoundingBoxFinder gefinder(bb, bb.maxlen()/3);
	for (auto& e: grid.vedges){
		gefinder.addentry(BoundingBox(*e->pfirst(), *e->plast()));
	}
	for (auto& fe: force_edges){
		for (auto isus: gefinder.suspects(BoundingBox(*fe.first, *fe.second))){
			auto& e = grid.vedges[isus];
			double ksieta[2];
			if (e->pfirst()->id >=0 || e->plast()->id >=0) continue;
			if (SectCross(*fe.first, *fe.second, *e->pfirst(), *e->plast(), ksieta)){
				bad_edges.push_back(isus);
			}
			
		}
	}
	Grid::Algos::RemoveEdges(grid, bad_edges);
	
	//get cell index for each split pair
	//and local indicies
	aa::constant_ids_pvec(grid.vvert, -1);
	for (int i=0; i<force_edges.size(); ++i){
		force_edges[i].first->id = i;
		force_edges[i].second->id = i;
	}
	vector<int> icells, lnode1, lnode2;
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
	Grid::Algos::CutCellDims(grid, 4);
}
};

GridData Mesher::TUnstructuredTriangleRecomb::_run(const Contour::Tree& source, const CoordinateMap2D<double>& embedded){
	GridData ret = gmsh_builder(source, embedded, 1, callback);

	callback->step_after(10, "Recombination check");
	//as a result of recombination some narrow reversed boundary triangles may occur.
	//here we try to fix it by merging with adjacent inner cells.
	recomb_heal(ret);
	//now we should check constraints lay on grid edges
	//since this feature can be not satisfied after recombination.
	EdgeData ec;
	for (auto n: source.nodes) if (n->isdetached()){
		for (auto e: n->contour){
			if (source.whereis(e->center()) == INSIDE) ec.push_back(e);
		}
	}
	guarantee_edges(ret, ec);
	
	return ret;
}

GridData Mesher::TUnstructuredTriangleRecomb::_run(const Contour::Tree& source){
	return _run(source, CoordinateMap2D<double>());
}
