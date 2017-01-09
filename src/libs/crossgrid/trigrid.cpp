#include "Gmsh.h"
#include "GModel.h"
#include "MVertex.h"
#include "MElement.h"
#include "addalgo.hpp"
#include "trigrid.h"
#include "nan_handler.h"
#include "procgrid.h"
#include "algos.hpp"
#include "constructor.hpp"
#include "treverter2d.hpp"

void TriGrid::FillFromGModel(void* gmod){
	GModel* m = static_cast<GModel*>(gmod);
	
	//--- build mesh
	//!! gmsh 2.11 doesn't work correctly without this line
	//   if chararcteristic mesh size is not 1.0
	auto bb = m->bounds();
	GmshSetBoundingBox(bb.min()[0], bb.max()[0], bb.min()[1], bb.max()[1], 0, 0);
	//GmshWriteFile("FF.opt");
	//m->writeGEO("gmsh_geo.geo");
	//shot down nan checks because gmsh has 1/0 operations in postprocessing
	NanSignalHandler::StopCheck();
	m->mesh(2);
	NanSignalHandler::StartCheck();
	//m->writeMSH("gmsh_msh.msh");
	//m->writeVTK("gmsh_msh.vtk");
	
	//--- extract mesh from gmsh
	std::map<MVertex*, GridPoint*> vrt;
	for (auto fit = m->firstFace(); fit!=m->lastFace(); ++fit){ 
		GFace* f = *fit;
		for (int i=0; i<f->getNumMeshElements(); ++i){
			auto e = f->getMeshElement(i);
			std::vector<MVertex*> vs;
			e->getVertices(vs);
			auto newcell = aa::add_shared(cells, Cell());
			for (auto v: vs){ 
				auto fnd = vrt.find(v); 
				if (fnd==vrt.end()){
					auto newp = aa::add_shared(points, GridPoint(v->x(), v->y()));
					auto ins = vrt.emplace(v, newp);
					fnd = ins.first;
				}
				add_point_to_cell(newcell, fnd->second);
			}
		}
	}

	//indicies
	set_indicies();
	//force positive triangles (necessary)
	force_cells_ordering();
}

TriGrid::TriGrid(const ContoursCollection& cont, const vector<double>& lc, double density){
	//build mesh in gmsh
	//2 - auto, 5 - delaunay, 6 - frontal
	GmshSetOption("Mesh", "Algorithm", 2.0);

	//very rough density parameter implementation: doesn't work in gmsh 2.9.1
	//double maxlc = 0;
	//for (auto v: lc) if (v>maxlc) maxlc=v;
	//GmshSetOption("Mesh", "CharacteristicLengthMax", (1.0-0.49*density)*maxlc); 

	GModel m;
	m.setFactory("Gmsh");
	//add points
	auto lcit = lc.begin();
	std::map<const Point*, GVertex*> verticies;
	for (auto i=0; i<cont.n_cont(); ++i){
		auto c = cont.get_contour(i);
		for (int j=0; j<c->n_points(); ++j){
			auto p = c->get_point(j);
			verticies[p] = m.addVertex(p->x, p->y, 0, *lcit++);
		}
	}
	//add edges and faces: each gmsh face should be single connected
	auto add_contour_edges = [&m, &verticies](std::vector<GEdge*>& e, const PContour* c){
		for (int j=0; j<c->n_points(); ++j){
			auto p = c->get_point(j);
			auto pn= c->get_point(j+1);
			e.push_back(m.addLine(verticies[p], verticies[pn]));
		}
	};
	std::vector<GFace*> fc;
	for (int i=0; i<cont.n_cont(); ++i) if (cont.is_inner(i)){
		std::vector<GEdge*> eds;
		auto c = cont.get_contour(i);
		add_contour_edges(eds, c);
		for (auto& oc: cont.get_childs(i)){
			add_contour_edges(eds, oc);
		}
		fc.push_back(m.addPlanarFace({eds}));
	}

	FillFromGModel(&m);
}

vector<HM2D::Contour::Tree>
TriGrid::ConstraintsPreproc(const HM2D::Contour::Tree& cont, 
		const ShpVector<HM2D::EdgeData>& constraints){
	vector<HM2D::Contour::Tree> ret;
	//1) sort out all inner contours
	for (auto rc: cont.nodes) if (rc->isbound()){
		if (rc->level % 2 == 0){
			ret.push_back(HM2D::Contour::Tree());
			auto& et = ret.back();
			et.add_contour(rc->contour);
			for (auto cc: rc->children){
				et.add_contour(cc.lock()->contour);
			}
		}
	}
	//2) check where constraint lies and add it to
	//   one of extended trees
	for (auto& c: constraints){
		for (auto p: HM2D::Contour::OrderedPoints(*c)){
			bool found = false;
			for (auto& et: ret){
				if (et.whereis(*p) != OUTSIDE){
					et.add_detached_contour(*c);
					found = true;
					break;
				}
			}
			if (found) break;
		}
	}

	/*

	//3) if point of open contour lies on bounding contours
	//   -> add this point to bounding contour
	auto split_edge = [](HM2D::EdgeData* cont, int ind, Point* point){
		Point* p1 = cont->edge(ind)->pstart;
		Point* p2 = cont->edge(ind)->pend;
		//equal nodes should have equal address
		if (*p1 == *point || *p2 == *point) {
			Point* p = (*p1 == *point) ? p1 : p2;
			auto f = cont->pinfo(p);
			if (f.eprev && f.eprev->pstart == p) f.eprev->pstart = point;
			if (f.eprev && f.eprev->pend == p) f.eprev->pend = point;
			if (f.enext && f.enext->pstart == p) f.enext->pstart = point;
			if (f.enext && f.enext->pend == p) f.enext->pend = point;
			return;
		}
		bool dircorrect = cont->correctly_directed_edge(ind);
		cont->RemoveAt({ind});
		auto e1 = std::make_shared<HMCont2D::Edge>(p1, point);
		auto e2 = std::make_shared<HMCont2D::Edge>(point, p2);
		if (dircorrect) cont->AddAt(ind, {e1, e2});
		else cont->AddAt(ind, {e2, e1});
	};
	for (auto& et: ret){
		for (auto& oc: et.open_contours){
			for (auto p: oc->all_points()){
				for (auto& bc: et.nodes){
					auto ca = bc->coord_at(*p);
					if (ISZERO(std::get<4>(ca))){
						split_edge(bc.get(), std::get<2>(ca), p);
						break;
					}
				}
			}
		}
	}
	*/

	return ret;
}

TriGrid::TriGrid(const HM2D::Contour::Tree& cont, 
		const ShpVector<HM2D::EdgeData>& constraints,
		double h){
	FillFromTree(cont, constraints, {}, std::map<Point*,double>(), h);
}

TriGrid::TriGrid(const HM2D::Contour::Tree& cont, 
		const ShpVector<HM2D::EdgeData>& constraints,
		const std::map<Point*, double>& w, double h){
	FillFromTree(cont, constraints, {}, w, h);
}

TriGrid::TriGrid(const HM2D::Contour::Tree& cont, 
		const ShpVector<HM2D::EdgeData>& constraints,
		const std::vector<double>& emb_points){
	std::map<Point*, double> w;

	vector<Point> ep; ep.reserve(emb_points.size()/3);
	for (int i=0; i<emb_points.size()/3; ++i){
		ep.push_back(Point(emb_points[3*i], emb_points[3*i+1]));
		w[&ep.back()] = emb_points[3*i+2];
	}

	FillFromTree(cont, constraints, ep, w, 0);
}

void TriGrid::CrossesProcessing(
		HM2D::Contour::Tree& cont, 
		ShpVector<HM2D::EdgeData>& constraints,
		std::map<Point*, double>& w,
		double h){
	auto getpw = [&](Point* p)->double{
		auto fnd = w.find(p);
		return (fnd==w.end())?h:fnd->second;
	};
	std::set<Point*> cross_points;
	auto treat_conts = [&](HM2D::EdgeData& c1, HM2D::EdgeData& c2){
		auto crosses = HM2D::Contour::Algos::CrossAll(c1, c2);
		for (auto& c: crosses){
			auto res1 = HM2D::Contour::GuaranteePoint(c1, std::get<1>(c));
			auto res2 = HM2D::Contour::GuaranteePoint(c2, std::get<1>(c));
			auto p1 = std::get<1>(res1);
			auto p2 = std::get<1>(res2);
			if (p1 == p2) continue;
			//substitute res2 with res1 pointers
			auto info1 = HM2D::Contour::PInfo(c1, p1.get());
			auto info2 = HM2D::Contour::PInfo(c2, p2.get());
			if (info2.eprev){
				if (info2.eprev->first() == p2) info2.eprev->vertices[0] = p1;
				else info2.eprev->vertices[1] = p1;
			}
			if (info2.enext){
				if (info2.enext->first() == p2) info2.enext->vertices[0] = p1;
				else info2.enext->vertices[1] = p1;
			}
			//set largest weight
			double cw = 0;
			if (!std::get<0>(res1)) cw = std::max(cw, getpw(p1.get()));
			else{
				if (info1.pprev != 0) cw = std::max(cw, getpw(info1.pprev.get()));
				if (info1.pnext != 0) cw = std::max(cw, getpw(info1.pnext.get()));
			}
			if (!std::get<0>(res2)) cw = std::max(cw, getpw(p2.get()));
			else{
				if (info2.pprev != 0) cw = std::max(cw, getpw(info2.pprev.get()));
				if (info2.pnext != 0) cw = std::max(cw, getpw(info2.pnext.get()));
			}
			w[p1.get()] = cw;
		}
	};
	//constraints vs constraints crosses
	for (int i=0; i<constraints.size(); ++i){
		for (int j=i+1; j<constraints.size(); ++j){
			auto c1 = constraints[i];
			auto c2 = constraints[j];
			treat_conts(*c1, *c2);
		}
	}
	//constraints vs contours
	for (auto& cns: constraints){
		for (auto& tree_cont: cont.nodes){
			treat_conts(tree_cont->contour, *cns);
		}
	}
}

void TriGrid::FillFromTree(
		const HM2D::Contour::Tree& cont_, 
		const ShpVector<HM2D::EdgeData>& constraints_,
		const vector<Point>& emb_points,
		const std::map<Point*, double>& w_,
		double h,
		bool recomb){
	if (cont_.nodes.size() == 0) return; 
	HM2D::Contour::R::RevertTree rt(cont_);

	//treat default size
	if (h<=0) h = 2*HM2D::BBox(cont_.alledges()).lendiag();

	//modify input data with respect to crosses with constraints
	HM2D::Contour::Tree cont = cont_;
	ShpVector<HM2D::EdgeData> constraints = constraints_;
	std::map<Point*, double> w = w_;
	if (constraints.size() > 0){
		cont = HM2D::Contour::Tree::DeepCopy(cont_);
		for (auto& c: constraints) {
			shared_ptr<HM2D::EdgeData> ne(new HM2D::EdgeData());
			DeepCopy(*c, *ne);
			std::swap(ne, c);
		}
		CrossesProcessing(cont, constraints, w, h);
	}

	//part tree by doubly connected ones and link each constraint with it.
	//returns ExtendedTree
	auto ap = ConstraintsPreproc(cont, constraints);

	//build mesh in gmsh
	//2 - auto, 5 - delaunay, 6 - frontal, 8 - delaunay for quads
	GmshSetOption("Mesh", "Algorithm", 2.0);
	GModel m;
	m.setFactory("Gmsh");

	//add points
	vector<const Point*> allpoints;
	for (auto& x: ap)
	for (auto& _t: HM2D::AllVertices(x.alledges())){
		allpoints.push_back(_t.get());
	}

	std::map<const Point*, GVertex*> verticies;
	for (auto& p: allpoints){
		auto wfnd = w.find(const_cast<Point*>(p));
		double hh = (wfnd == w.end()) ? h : wfnd->second;
		verticies[p] = m.addVertex(p->x, p->y, 0, hh);
	}
	

	//add edges and faces: each gmsh face should be single connected
	auto add_contour_edges = [&m, &verticies](const HM2D::EdgeData& c,
			vector<GEdge*>& e){
		HM2D::VertexData op = HM2D::Contour::OrderedPoints(c);
		for (int i=0; i<op.size()-1; ++i){
			auto p = op[i].get();
			auto pn = op[i+1].get();
			e.push_back(m.addLine(verticies[p], verticies[pn]));
		}
	};

	NanSignalHandler::StopCheck();
	//add edges
	std::vector<GFace*> fc;
	for (auto& ec: ap){
		std::vector<GEdge*> eds;
		//inner contour
		add_contour_edges(ec.roots()[0]->contour, eds);
		//outer contours
		for (auto& child: ec.roots()[0]->children){
			add_contour_edges(child.lock()->contour, eds);
		}
		//assemble face
		fc.push_back(m.addPlanarFace({eds}));
		//constraints
		eds.clear();
		for (auto c: ec.detached_contours()){
			add_contour_edges(c->contour, eds);
		}
		for (auto e: eds) fc.back()->addEmbeddedEdge(e);
	}
	NanSignalHandler::StartCheck();

	//add embedded points
	for (auto& p: emb_points){
		auto wfnd = w.find(const_cast<Point*>(&p));
		double hh = (wfnd == w.end()) ? h : wfnd->second;
		auto added = m.addVertex(p.x, p.y, 0, hh);
		//find face containing point
		for (int i=0; i<ap.size(); ++i){
			if (ap[i].whereis(p) == INSIDE) fc[i]->addEmbeddedVertex(added);
			break;
		}
	}

	if (recomb){
		//build 1d mesh explicitly without recombination because
		//otherwise gmsh make boundaries twice as fine
		m.mesh(1);
		//usage of delaunay for quads gives worse results for non-regular areas
		//hence using auto algorithm
		GmshSetOption("Mesh", "Algorithm", 2.0);
		GmshSetOption("Mesh", "RecombinationAlgorithm", 1.0);
		for (auto& f: fc) f->meshAttributes.recombine = 1.0;
	}

	FillFromGModel(&m);

	if (recomb){
		// if all nodes of quad grid are still trianlge than most likely builder
		// has failed. So we use another algorithm
		bool has4=false;
		for (auto c: cells) if (c->dim() == 4) {has4=true; break; }
		if (!has4){
			clear();
			GmshSetOption("Mesh", "Algorithm", 2.0);
			GmshSetOption("Mesh", "RecombinationAlgorithm", 0.0);
			for (auto& f: fc) f->meshAttributes.recombine = 1.0;
			FillFromGModel(&m);
		}
		//as a result of recombination some narrow reversed boundary triangles may occur.
		//here we try to fix it by merging with adjacent inner cells.
		recomb_heal();
		GGeom::Repair::CellsTo34(*this);
		//now we should check constraints lay on grid edges
		//since this feature can be not satisfied after recombination.
		vector<HM2D::Edge*> alledges;
		for (auto& ev: constraints)
		for (auto& e: (*ev)) alledges.push_back(e.get());
		guarantee_edges(alledges);
	}
}

void TriGrid::recomb_heal(){
	auto process = [&](int ic1, int ic2, int n1, int n2){
		Cell* c1 = cells[ic1].get();
		Cell* c2 = cells[ic2].get();
		std::reverse(c2->points.begin(), c2->points.end());
		int loc1 = -1, loc2 = -1;
		for (int j=0; j<c1->dim(); ++j){
			if (c1->get_point(j)->get_ind() == n1){
				loc1 = j; break;
			}
		}
		for (int j=0; j<c2->dim(); ++j){
			if (c2->get_point(j)->get_ind() == n1){
				loc2 = j; break;
			}
		}
		assert(loc1 >= 0 && loc2 >= 0);
		assert(c1->get_point(loc1+1)->get_ind() == n2);
		std::rotate(c2->points.begin(), c2->points.begin()+loc2, c2->points.end());
		c1->points.insert(c1->points.begin()+loc1+1, c2->points.begin()+1, c2->points.end()-1);
		GGeom::Modify::RemoveCells(*this, {c2});
		return !c1->has_self_crosses();
	};

	for (int tries=0; tries<100; ++tries){
		std::map<std::pair<int, int>, int> used_edges;
		set_indicies();
		for (int ic=0; ic<cells.size(); ++ic){
			Cell* c=cells[ic].get();
			for (int i=0; i<c->dim(); ++i){
				int p1 = c->get_point(i)->get_ind();
				int p2 = c->get_point(i+1)->get_ind();
				auto er = used_edges.emplace(std::make_pair(p1, p2), ic);
				if (er.second==false){
					if (process(ic, er.first->second, p1, p2)) goto NEXT_TRY;
					else goto ERR_OUT;
				}
			}
		}
		//no changes in current try => grid is fine.
		return;
NEXT_TRY:
		continue; 
	}
ERR_OUT:
	throw std::runtime_error("cannot restore correct grid from gmsh output");
}

void TriGrid::guarantee_edges(const vector<HM2D::Edge*>& ed){
	auto cf = GGeom::Info::CellFinder(this, 30, 30);
	std::map<const Cell*, int> divide_cells;
	vector<Point> midpoints; midpoints.reserve(ed.size());
	for (auto& ev: ed) midpoints.push_back(ev->center());
	for (size_t i=0; i<midpoints.size(); ++i){
		auto candcells = cf.CellCandidates(midpoints[i]);
		Point* pstart = ed[i]->first().get();
		Point* pend = ed[i]->last().get();
		for (auto cand: candcells) if (cand->dim() == 4){
			auto cc = GGeom::Info::CellContour(*this, cand->get_ind());
			if (HM2D::Contour::WhereIs(cc, midpoints[i]) == INSIDE){
				int ep1=-1, ep2=-1;
				for (int j=0; j<cand->dim(); ++j){
					auto p1 = cand->get_point(j);
					if (*p1 == *pstart) ep1 = j;
					if (*p1 == *pend) ep2 = j;
				}
				if (ep1>=0 && ep2>=0){
					if (ep2 < ep1) std::swap(ep1, ep2);
					if ((ep2-ep1) % 2 == 0){
						auto fnd = divide_cells.find(cand);
						if (fnd != divide_cells.end()){
							fnd->second = -1;
						} else {
							divide_cells.emplace(cand, ep1);
						}
					}
				}
				break;
			}
		}
	}
	for (auto kv: divide_cells) if (kv.second>=0){
		GridPoint* p0 = kv.first->points[(kv.second) % 4];
		GridPoint* p1 = kv.first->points[(kv.second+1) % 4];
		GridPoint* p2 = kv.first->points[(kv.second+2) % 4];
		GridPoint* p3 = kv.first->points[(kv.second+3) % 4];
		cells[kv.first->get_ind()]->points = {p0, p1, p2};
		auto newcell = aa::add_shared(cells, Cell());
			newcell->points = {p0, p2, p3};
	}
	set_cell_indicies();
}

shared_ptr<TriGrid> TriGrid::FromGmshGeo(const char* fn){
	GmshSetOption("Mesh", "Algorithm", 2.0);
	GModel m;
	m.setFactory("Gmsh");
	m.load(fn);
	shared_ptr<TriGrid> ret(new TriGrid());
	ret->FillFromGModel(&m);
	return ret;
}

std::set<Edge>& TriGrid::edges() const{
	if (_edges.size()==0) _edges = get_edges();
	return _edges;
}
std::map<GridPoint*, vector<GridPoint*>>& TriGrid::nodenodeI() const{
	if (_nodenodeI.size()==0){
		//1) get nodenode
		auto& eds = edges();
		for (auto& e: eds){
			auto it1 = _nodenodeI.emplace(points[e.p1].get(), vector<GridPoint*>());
			it1.first->second.push_back(points[e.p2].get());
			auto it2 = _nodenodeI.emplace(points[e.p2].get(), vector<GridPoint*>());
			it2.first->second.push_back(points[e.p1].get());
		}
		//2) remove boundary
		for (auto& e: eds) if (e.is_boundary()){
			auto it1 = _nodenodeI.find(points[e.p1].get());
			auto it2 = _nodenodeI.find(points[e.p2].get());
			if (it1!=_nodenodeI.end()) _nodenodeI.erase(it1);
			if (it2!=_nodenodeI.end()) _nodenodeI.erase(it2);
		}
	}
	return _nodenodeI;
}

ShpVector<Point> TriGrid::ref_points(const vector<double>& dists, double density) const{
	auto ret=ShpVector<Point>();
	for (auto e: edges()){
		if (!e.is_boundary()){
			auto p1 = get_point(e.p1), p2 = get_point(e.p2);
			double len = Point::dist(*p1, *p2);
			auto ksi = RefineSection(dists[e.p1], dists[e.p2], len, density);
			for (auto k: ksi) aa::add_shared(ret, Point::Weigh(*p1, *p2, k/len)); 
		}
	}
	return ret;
}

void TriGrid::smooth(double w){
	for (auto& nn: nodenodeI()){
		Point wav(0,0);
		for (auto& an: nn.second) wav+=*an;
		wav*=(w/nn.second.size());
		(*nn.first)*=(1-w); (*nn.first)+=wav;
	}
}

vector<Point> TriGrid::cell_centers() const{
	vector<Point> ret;
	for (auto& c: cells){
		auto p = Point(0,0);
		p += *c->get_point(0);
		p += *c->get_point(1);
		p += *c->get_point(2);
		p/=3;
		ret.push_back(p);
	}
	return ret;
}

vector<double> TriGrid::cell_areas() const{
	vector<double> ret;
	for (auto& c: cells){
		ret.push_back(triarea(*c->get_point(0), *c->get_point(1), *c->get_point(2)));
	}
	return ret;
}

shared_ptr<TriGrid>
TriGrid::TriangulateArea(const vector<Point>& pts, double h){
	Contour c(pts);
	ContoursCollection cc({c});
	vector<double> lc (c.n_points(), h);
	shared_ptr<TriGrid> ret(new TriGrid(cc, lc, 7.0));
	return ret;
}

shared_ptr<TriGrid>
TriGrid::TriangulateArea(const vector<vector<Point>>& pts, double h){
	vector<PContour> c;
	for (auto& x: pts){
		PContour pc;
		for (auto& x2: x) pc.add_point(const_cast<Point*>(&x2));
		c.push_back(pc);
	}
	ContoursCollection cc(c);
	int np = 0;
	for (auto& x: pts) np += x.size();
	vector<double> lc (np, h);
	shared_ptr<TriGrid> ret(new TriGrid(cc, lc, 7.0));
	return ret;
}

shared_ptr<TriGrid>
TriGrid::TriangulateAreaConstrained(const vector<vector<Point>>& bnd,
		const vector<vector<Point>>& cns, double h){
	ShpVector<HM2D::EdgeData> contours;
	ShpVector<HM2D::EdgeData> constraints;
	for (auto& p: bnd) aa::add_shared(
			contours,
			HM2D::Contour::Constructor::FromPoints(p, true)
	);

	for (auto& p: cns) aa::add_shared(
			constraints,
			HM2D::Contour::Constructor::FromPoints(p, false)
	);

	ShpVector<HM2D::EdgeData> ccontours;
	for (auto v: contours) ccontours.push_back(v);
	ShpVector<HM2D::EdgeData> cconstraints;
	for (auto v: constraints) cconstraints.push_back(v);


	HM2D::Contour::Tree tree;
	for (auto& p: ccontours){
		tree.add_contour(*p);
	}


	//build
	shared_ptr<TriGrid> ret(new TriGrid(tree, cconstraints, h));
	return ret;
}

shared_ptr<TriGrid>
TriGrid::TriangulateArea(const HM2D::Contour::Tree& cont, double h){
	return shared_ptr<TriGrid>(new TriGrid(cont, {}, h));
}

shared_ptr<TriGrid>
TriGrid::TriangulateArea(const HM2D::Contour::Tree& cont, const std::map<Point*, double>& w, double h){
	return shared_ptr<TriGrid>(new TriGrid(cont, {}, w, h));
}


GridGeom QuadGrid(const HM2D::Contour::Tree& cont, 
		const ShpVector<HM2D::EdgeData>& constraints,
		const std::vector<double>& emb_points){
	std::map<Point*, double> w;

	vector<Point> ep; ep.reserve(emb_points.size()/3);
	for (int i=0; i<emb_points.size()/3; ++i){
		ep.push_back(Point(emb_points[3*i], emb_points[3*i+1]));
		w[&ep.back()] = emb_points[3*i+2];
	}

	TriGrid g;
	g.FillFromTree(cont, constraints, ep, w, 0, true);
	return GridGeom(std::move(g));
}

shared_ptr<GridGeom> QuadrangulateArea(const HM2D::Contour::Tree& cont,
		const std::map<Point*, double>& w, double h){
	TriGrid g;
	g.FillFromTree(cont, {}, {}, w, 2, true);
	shared_ptr<GridGeom> ret(new GridGeom(std::move(g)));
	return ret;
}
