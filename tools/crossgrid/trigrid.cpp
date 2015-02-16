#include "addalgo.hpp"
#include "trigrid.h"
#include "gmsh/Gmsh.h"
#include "gmsh/GModel.h"
#include "gmsh/MVertex.h"
#include "gmsh/MElement.h"

TriGrid::TriGrid(const vector<PContour>& cont, const vector<double>& lc){
	//build mesh in gmsh
	GmshInitialize();
	GModel m;
	m.setFactory("Gmsh");
	auto lcit = lc.begin();
	std::vector<GEdge*> edges;
	for (auto c: cont){
		std::vector<GVertex*> verticies;
		for (int i=0; i<c.n_points(); ++i){
			auto p = c.get_point(i);
			verticies.push_back(m.addVertex(p->x, p->y, 0, *lcit++));
		}
		auto pprev = verticies.back();
		for (auto p: verticies){
			edges.push_back(m.addLine(pprev, p));
			pprev = p;
		}
	}
	std::vector<std::vector<GEdge*>> loop;
	loop.push_back(edges);
	GFace *f = m.addPlanarFace(loop);
	m.mesh(2);
	//--- extract mesh from gmsh
	std::map<MVertex*, GridPoint*> vrt;
	//cells
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

	GmshFinalize();

	//indicies
	set_indicies();
	//force positive triangles (necessary)
	force_cells_ordering();
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

shp_vector<Point> TriGrid::ref_points(const vector<double>& dists, double density) const{
	auto ret=shp_vector<Point>();
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

