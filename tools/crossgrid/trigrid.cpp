#include "gmsh/Gmsh.h"
#include "gmsh/GModel.h"
#include "gmsh/MVertex.h"
#include "gmsh/MElement.h"
#include "addalgo.hpp"
#include "trigrid.h"

TriGrid::TriGrid(const ContoursCollection& cont, const vector<double>& lc, double density){
	//build mesh in gmsh
	//2 - auto, 5 - delauney, 6 - frontal
	GmshSetOption("Mesh", "Algorithm", 6.0);

	//very rough density parameter implementation
	double maxlc = 0;
	for (auto v: lc) if (v>maxlc) maxlc=v;
	GmshSetOption("Mesh", "CharacteristicLengthMax", (1.0-0.49*density)*maxlc); 

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

	//--- build mesh
	//m.writeGEO("gmsh_geo.geo");
	m.mesh(2);
	//m.writeMSH("gmsh_msh.msh");
	
	//--- extract mesh from gmsh
	std::map<MVertex*, GridPoint*> vrt;
	for (auto& f: fc){
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

