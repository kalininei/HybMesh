#include "Gmsh.h"
#include "GModel.h"
#include "MVertex.h"
#include "MElement.h"
#include "addalgo.hpp"
#include "trigrid.h"
#include "fileproc.h"

void TriGrid::FillFromGModel(void* gmod){
	GModel* m = static_cast<GModel*>(gmod);
	
	//--- build mesh
	//m->writeGEO("gmsh_geo.geo");
	m->mesh(2);
	//m->writeMSH("gmsh_msh.msh");
	
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
	//##############################################
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
		Point p = *c->first();
		for (auto& et: ret){
			if (!et.IsWithout(p)){
				shared_ptr<HMCont2D::Contour> sc(c);
				et.AddOpenContour(sc);
				break;
			}
		}
	}

	//3) if point of open contour lie on bounding contours
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

	return ret;
}

TriGrid::TriGrid(const HMCont2D::ContourTree& cont, 
		const ShpVector<HMCont2D::Contour>& constraints,
		double h){
	FillFromTree(cont, constraints, std::map<Point*,double>(), h);
}

TriGrid::TriGrid(const HMCont2D::ContourTree& cont, 
		const ShpVector<HMCont2D::Contour>& constraints,
		const std::map<Point*, double>& w, double h){
	FillFromTree(cont, constraints, w, h);
}

void TriGrid::FillFromTree(const HMCont2D::ContourTree& cont, 
		const ShpVector<HMCont2D::Contour>& constraints,
		const std::map<Point*, double>& w, double h){
	auto ap = ConstraintsPreproc(cont, constraints);

	//build mesh in gmsh
	//2 - auto, 5 - delaunay, 6 - frontal
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

	FillFromGModel(&m);
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
	for (auto c: cells){
		auto p = Point(0,0);
		p += *c->get_point(0);
		p += *c->get_point(1);
		p += *c->get_point(2);
		p/=3;
		ret.push_back(p);
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


