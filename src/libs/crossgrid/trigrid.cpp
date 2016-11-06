#include "Gmsh.h"
#include "GModel.h"
#include "MVertex.h"
#include "MElement.h"
#include "addalgo.hpp"
#include "trigrid.h"
#include "nan_handler.h"
#include "procgrid.h"

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

vector<HMCont2D::ExtendedTree>
TriGrid::ConstraintsPreproc(const HMCont2D::ContourTree& cont, 
		const ShpVector<HMCont2D::Contour>& constraints){
	vector<HMCont2D::ExtendedTree> ret;
	//1) sort out all inner contours
	for (auto rc: cont.nodes){
		if (HMCont2D::Contour::Area(*rc) > 0){
			ret.push_back(HMCont2D::ExtendedTree());
			auto& et = ret.back();
			shared_ptr<HMCont2D::Contour> rc2(new HMCont2D::Contour(*rc));
			et.AddContour(rc2);
			for (auto cc: rc->children){
				shared_ptr<HMCont2D::Contour> cc2(new HMCont2D::Contour(*cc));
				et.AddContour(cc2);
			}
		}
	}
	//2) check where constraint lies and add it to
	//   one of extended trees
	for (auto& c: constraints){
		for (auto p: c->ordered_points()){
			bool found = false;
			for (auto& et: ret){
				if (!et.IsWithout(*p)){
					shared_ptr<HMCont2D::Contour> sc(c);
					et.AddOpenContour(sc);
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
	auto split_edge = [](HMCont2D::Contour* cont, int ind, Point* point){
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

TriGrid::TriGrid(const HMCont2D::ContourTree& cont, 
		const ShpVector<HMCont2D::Contour>& constraints,
		double h){
	FillFromTree(cont, constraints, {}, std::map<Point*,double>(), h);
}

TriGrid::TriGrid(const HMCont2D::ContourTree& cont, 
		const ShpVector<HMCont2D::Contour>& constraints,
		const std::map<Point*, double>& w, double h){
	FillFromTree(cont, constraints, {}, w, h);
}

TriGrid::TriGrid(const HMCont2D::ContourTree& cont, 
		const ShpVector<HMCont2D::Contour>& constraints,
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
		HMCont2D::ContourTree& cont, 
		ShpVector<HMCont2D::Contour>& constraints,
		std::map<Point*, double>& w,
		HMCont2D::PCollection& apnt,
		double h){
	auto getpw = [&](Point* p)->double{
		auto fnd = w.find(p);
		return (fnd==w.end())?h:fnd->second;
	};
	std::set<Point*> cross_points;
	auto treat_conts = [&](HMCont2D::Contour& c1, HMCont2D::Contour& c2){
		auto crosses = HMCont2D::Algos::CrossAll(c1, c2);
		for (auto& c: crosses){
			auto res1 = c1.GuaranteePoint(std::get<1>(c), apnt);
			auto res2 = c2.GuaranteePoint(std::get<1>(c), apnt);
			auto p1 = std::get<1>(res1);
			auto p2 = std::get<1>(res2);
			if (p1 == p2) continue;
			//substitute res2 with res1 pointers
			auto info1 = c1.pinfo(p1);
			auto info2 = c2.pinfo(p2);
			if (info2.eprev){
				if (info2.eprev->pstart == p2) info2.eprev->pstart = p1;
				else info2.eprev->pend = p1;
			}
			if (info2.enext){
				if (info2.enext->pstart == p2) info2.enext->pstart = p1;
				else info2.enext->pend = p1;
			}
			//set largest weight
			double cw = 0;
			if (!std::get<0>(res1)) cw = std::max(cw, getpw(p1));
			else{
				if (info1.pprev != 0) cw = std::max(cw, getpw(info1.pprev));
				if (info1.pnext != 0) cw = std::max(cw, getpw(info1.pnext));
			}
			if (!std::get<0>(res2)) cw = std::max(cw, getpw(p2));
			else{
				if (info2.pprev != 0) cw = std::max(cw, getpw(info2.pprev));
				if (info2.pnext != 0) cw = std::max(cw, getpw(info2.pnext));
			}
			w[p1] = cw;
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
			treat_conts(*tree_cont, *cns);
		}
	}
}

void TriGrid::FillFromTree(
		const HMCont2D::ContourTree& cont_, 
		const ShpVector<HMCont2D::Contour>& constraints_,
		const vector<Point>& emb_points,
		const std::map<Point*, double>& w_,
		double h,
		bool recomb){
	//treat default size
	if (h<=0) h = 2*HMCont2D::ECollection::BBox(cont_).lendiag();

	//modify input data with respect to crosses with constraints
	HMCont2D::ContourTree cont = cont_;
	ShpVector<HMCont2D::Contour> constraints = constraints_;
	std::map<Point*, double> w = w_;
	HMCont2D::PCollection apnt;
	if (constraints.size() > 0){
		for (auto& c: cont.nodes) c->Reallocate();
		cont.ReloadEdges();
		for (auto& c: constraints) c->Reallocate();
		CrossesProcessing(cont, constraints, w, apnt, h);
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
	for (auto& x: ap){
		auto _t = x.all_points();
		std::copy(_t.begin(), _t.end(), std::back_inserter(allpoints));
	}

	std::map<const Point*, GVertex*> verticies;
	for (auto& p: allpoints){
		auto wfnd = w.find(const_cast<Point*>(p));
		double hh = (wfnd == w.end()) ? h : wfnd->second;
		verticies[p] = m.addVertex(p->x, p->y, 0, hh);
	}
	

	//add edges and faces: each gmsh face should be single connected
	auto add_contour_edges = [&m, &verticies](const HMCont2D::Contour& c,
			vector<GEdge*>& e){
		vector<Point*> op = c.ordered_points();
		for (int i=0; i<op.size()-1; ++i){
			const Point* p = op[i];
			const Point* pn = op[i+1];
			e.push_back(m.addLine(verticies[p], verticies[pn]));
		}
	};

	//add edges
	std::vector<GFace*> fc;
	for (auto& ec: ap){
		std::vector<GEdge*> eds;
		//inner contour
		add_contour_edges(*ec.roots()[0], eds);
		//outer contours
		for (auto& child: ec.roots()[0]->children){
			add_contour_edges(*child, eds);
		}
		//assemble face
		fc.push_back(m.addPlanarFace({eds}));
		//constraints
		eds.clear();
		for (int i=ec.nodes.size(); i<ec.cont_count(); ++i){
			add_contour_edges(*ec.get_contour(i), eds);
		}
		for (auto e: eds) fc.back()->addEmbeddedEdge(e);
	}

	//add embedded points
	for (auto& p: emb_points){
		auto wfnd = w.find(const_cast<Point*>(&p));
		double hh = (wfnd == w.end()) ? h : wfnd->second;
		auto added = m.addVertex(p.x, p.y, 0, hh);
		//find face containing point
		for (int i=0; i<ap.size(); ++i){
			if (ap[i].IsWithin(p)) fc[i]->addEmbeddedVertex(added);
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

	// if all nodes of quad grid are still trianlge than most likely builder
	// has failed. So we use another algorithm
	if (recomb){
		bool has4=false;
		for (auto c: cells) if (c->dim() == 4) {has4=true; break; }
		if (!has4){
			clear();
			GmshSetOption("Mesh", "Algorithm", 2.0);
			GmshSetOption("Mesh", "RecombinationAlgorithm", 0.0);
			for (auto& f: fc) f->meshAttributes.recombine = 1.0;
			FillFromGModel(&m);
		}
	}
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
	ShpVector<HMCont2D::Container<HMCont2D::Contour>> contours;
	ShpVector<HMCont2D::Container<HMCont2D::Contour>> constraints;
	for (auto& p: bnd) aa::add_shared(
			contours,
			HMCont2D::Constructor::ContourFromPoints(p, true)
	);

	for (auto& p: cns) aa::add_shared(
			constraints,
			HMCont2D::Constructor::ContourFromPoints(p, false)
	);

	ShpVector<HMCont2D::Contour> ccontours;
	for (auto v: contours) ccontours.push_back(v);
	ShpVector<HMCont2D::Contour> cconstraints;
	for (auto v: constraints) cconstraints.push_back(v);


	HMCont2D::ContourTree tree;
	for (auto& p: ccontours){
		tree.AddContour(p);
	}


	//build
	shared_ptr<TriGrid> ret(new TriGrid(tree, cconstraints, h));
	return ret;
}

shared_ptr<TriGrid>
TriGrid::TriangulateArea(const HMCont2D::ContourTree& cont, double h){
	return shared_ptr<TriGrid>(new TriGrid(cont, {}, h));
}

shared_ptr<TriGrid>
TriGrid::TriangulateArea(const HMCont2D::ContourTree& cont, const std::map<Point*, double>& w, double h){
	return shared_ptr<TriGrid>(new TriGrid(cont, {}, w, h));
}


GridGeom QuadGrid(const HMCont2D::ContourTree& cont, 
		const ShpVector<HMCont2D::Contour>& constraints,
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

shared_ptr<GridGeom> QuadrangulateArea(const HMCont2D::ContourTree& cont,
		const std::map<Point*, double>& w, double h){
	TriGrid g;
	g.FillFromTree(cont, {}, {}, w, 2, true);
	shared_ptr<GridGeom> ret(new GridGeom(std::move(g)));
	return ret;
}

namespace{
struct PebiBndPoint{
	int ind;
	shared_ptr<GridPoint> prev, next;
	shared_ptr<GridPoint> hprev, hnext;
};

std::vector<PebiBndPoint> assemble_bnd(const TriGrid& g){
	std::vector<PebiBndPoint> ret;
	auto cont = GGeom::Info::Contour(g);
	auto gpoints = GGeom::Info::SharePoints(g);
	for (auto c: cont.nodes){
		auto op = c->ordered_points();
		for (int i=0; i<op.size()-1; ++i){
			int im = (i==0)?op.size()-2:i-1;
			int cur = static_cast<const GridPoint*>(op[i])->get_ind();
			int prev = static_cast<const GridPoint*>(op[im])->get_ind();
			int next = static_cast<const GridPoint*>(op[i+1])->get_ind();
			ret.push_back(PebiBndPoint {cur, gpoints[prev], gpoints[next]});
			ret.back().hnext.reset(new GridPoint((*gpoints[cur]+*gpoints[next])/2.0));
			if (i!=0) ret.back().hprev = ret.end()[-2].hnext;
			if (i==op.size()-2) ret[ret.size()-1-i].hprev = ret.back().hnext;
		}
	}
	return ret;
}

std::vector<std::vector<int>> ordered_points_cells(const TriGrid& g, std::vector<PebiBndPoint>& bnd){
	std::vector<std::vector<int>> ret(g.n_points());
	auto edges=g.get_edges();
	std::vector<std::list<Edge>> point_edges(g.n_points());
	for (auto& e: edges){
		point_edges[e.p1].push_back(e);
		point_edges[e.p2].push_back(e);
	}
	for (int i=0; i<g.n_points(); ++i){
		auto& pe = point_edges[i];
		Edge estart(-1, -1);
		//looking for boundary edge
		for (auto it = pe.begin(); it!=pe.end(); ++it){
			int right_cell = (it->p1 == i)?it->cell_right:it->cell_left;
			if (right_cell < 0){ estart = *it; pe.erase(it); break; }
		}
		//if no boundary edges start from first
		if (estart.p1 == -1) {estart = pe.front(); pe.pop_front();}
		Edge ecur = estart;
		while (1){
			int cur_cell = (ecur.p1 == i)?ecur.cell_left:ecur.cell_right;
			if (cur_cell<0) break;
			ret[i].push_back(cur_cell);
			//find cur_edge as edge of cur_cell which is not ecur but contains i;
			if (pe.size() == 0) break;
			else for (auto it=pe.begin(); it!=pe.end(); ++it){
				//find edge which has a link to cur_cell
				if (it->cell_left == cur_cell || it->cell_right == cur_cell){
					ecur = *it; pe.erase(it); break;
				}
			}
		}
	}

	return ret;
}

Point calc_pebi(Point p1, Point p2, Point p3){
	auto line2p = [](double x1, double y1, double x2, double y2){
		double xc = (x1+x2)/2.0, yc = (y1+y2)/2.0;
		return std::array<double, 3> { x2-x1, y2-y1, -(x2-x1)*xc-(y2-y1)*yc };
	};
	auto l1 = line2p(p1.x, p1.y, p2.x, p2.y);
	auto l2 = line2p(p2.x, p2.y, p3.x, p3.y);
	double A[] = {l1[0], l1[1], l2[0], l2[1]};
	double det = A[0]*A[3]-A[1]*A[2];
	if (ISZERO(det)) throw std::runtime_error("pebi error");
	double B[] = {A[3]/det, -A[1]/det, -A[2]/det, A[0]/det};
	return Point(-B[0]*l1[2]-B[1]*l2[2], -B[2]*l1[2]-B[3]*l2[2]);
}

ShpVector<GridPoint> build_pebi_pts(const TriGrid& g){
	ShpVector<GridPoint> ret;
	vector<vector<int>> cc = g.cell_cell();
	auto within = [&](Point p, int cn)->bool{
		const Cell* c=g.get_cell(cn);
		auto &x1 = c->points[0]->x, &x2 = c->points[1]->x, &x3 = c->points[2]->x;
		auto &y1 = c->points[0]->y, &y2 = c->points[1]->y, &y3 = c->points[2]->y;
		double j11 = x2 - x1, j21 = x3 - x1;
		double j12 = y2 - y1, j22 = y3 - y1;
		double modj = (j22*j11 - j21*j12);
		Point ksieta;
		ksieta.x = ( j22*(p.x - x1) - j21*(p.y - y1))/modj;
		ksieta.y = (-j12*(p.x - x1) + j11*(p.y - y1))/modj;
		return (ksieta.x >= 0 && ksieta.x <= 1.0 && ksieta.y>=0 && ksieta.y<=1.0-ksieta.x);
	};
	for (int i=0; i<g.n_cells(); ++i){
		Point p1 = *g.get_cell(i)->points[0];
		Point p2 = *g.get_cell(i)->points[1];
		Point p3 = *g.get_cell(i)->points[2];
		Point p = calc_pebi(p1, p2, p3);

		//bad point is point which lies far away outside parent triangle.
		//Usually even out of area. Hence we snap it to parant triangle edge.
		bool bad_point=true;
		//within itself
		if (within(p, i)) bad_point = false;
		//within neighbours
		if (bad_point) for (int j=0; j<cc[i].size(); ++j){
			if (within(p, cc[i][j])) { bad_point = false; break; }
		}
		//within neighbours of neighbours
		if (bad_point) for (int j=0; j<cc[i].size(); ++j){
			int in = cc[i][j];
			for (int k=0; k<cc[in].size(); ++k)
				if (within(p, cc[in][k])) { bad_point = false; break; }
			if (!bad_point) break;
		}

		if (bad_point){
			double d1 = Point::meas_line(p, p1, p2);
			double d2 = Point::meas_line(p, p1, p3);
			double d3 = Point::meas_line(p, p2, p3);
			if (d1 < d2 && d1 < d3) p = (p1 + p2)/2.0;
			else if (d2<d1 && d2<d3) p = (p1+p3)/2.0;
			else p = (p2+p3)/2.0;
		}
		aa::add_shared(ret, GridPoint(p));
	}
	return ret;
}


}
GridGeom TriGrid::ToPeBi() const{
	//calculate pebi points
	ShpVector<GridPoint> pebi_points = build_pebi_pts(*this);

	//boundary information
	std::vector<PebiBndPoint> bnd_left_right = assemble_bnd(*this);

	//ordered points->cells table
	std::vector<std::vector<int>> points_cells = ordered_points_cells(*this, bnd_left_right);

	//asseble points
	ShpVector<GridPoint> rp(pebi_points.begin(), pebi_points.end());
	for (int i=0; i<bnd_left_right.size(); ++i){
		rp.push_back(bnd_left_right[i].hnext);
	}
	//assemble cells
	ShpVector<Cell> rc;
	//internal cells
	for (int i=0; i<n_points(); ++i){
		Cell* c=aa::add_shared(rc, Cell());
		for (auto& x: points_cells[i]) c->points.push_back(pebi_points[x].get());
	}
	//boundary segments
	for (auto& b: bnd_left_right){
		Cell* c = rc[b.ind].get();
		c->points.insert(c->points.begin(), b.hnext.get());
		c->points.push_back(b.hprev.get());
		double ksi;
		if (!isOnSection(*points[b.ind], *b.hnext, *b.hprev, ksi)){
			c->points.push_back(aa::add_shared(rp, *points[b.ind]));
		}
	}
	
	auto ret = GGeom::Constructor::FromData(rp, rc);
	//post processing: collapse reversed edges which provoke сell self intersection if possible
	std::set<const GridPoint*> bpts = get_bnd_points();
	for (int i=0; i<rc.size(); ++i){
		const Cell* c = rc[i].get();
		for (int k=0; k<c->dim(); ++k){
			GridPoint* p0 = const_cast<GridPoint*>(c->get_point(k));
			GridPoint* p1 = const_cast<GridPoint*>(c->get_point(k+1));
			int w = LinePointWhereIs(*points[i], *p0, *p1);
			if (w != 2) continue;
			if (bpts.find(p0) != bpts.end() || bpts.find(p1) != bpts.end()) continue;
			GridPoint* prev = const_cast<GridPoint*>(c->get_point(c->dim()+k-1));
			GridPoint* next = const_cast<GridPoint*>(c->get_point(k+2));
			if (*prev == *p0 || *next == *p1) continue;
			double ksieta[2];
			if (SectCross(*prev, *p0, *p1, *next, ksieta)){
				p0->set((*p0 + *p1)/2.0);
				p1->set(*p0);
			}
		}
	}
	//remove short edges
	GGeom::Repair::RemoveShortEdges(ret, 0.1);
	return ret;
}

