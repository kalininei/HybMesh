#include "vtkSmartPointer.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPolygon.h"
#include "vtkCellArray.h"
#include "vtkDelaunay2D.h"
#include "addalgo.hpp"
#include "trigrid.h"

vtkSmartPointer<vtkPoints> build_vtk_points(const vector<const Point*>& pts){
	vtkSmartPointer<vtkPoints> ret = vtkSmartPointer<vtkPoints>::New();
	for (auto p: pts) ret->InsertNextPoint(p->x, p->y, 0);
	return ret;
}

//perform constrained Delaunay 
std::vector<int> perform_triangulation(const vector<PContour>& cont, const shp_vector<Point>& pts){
	//1) copy points to vtk
	vector<const Point*> _pc;
	for (auto c: cont)
		for (int i=0; i<c.n_points(); ++i) _pc.push_back(c.get_point(i));
	for (auto p: pts) _pc.push_back(p.get());
	auto vtkPts = build_vtk_points(_pc);
	
	//2) build a poly data from points
	vtkSmartPointer<vtkPolyData> vtkPoly = vtkSmartPointer<vtkPolyData>::New();
	vtkPoly->SetPoints(vtkPts);
	
	//3) build a boundary cell array for each cont
	vtkSmartPointer<vtkCellArray> vtkBndCell = vtkSmartPointer<vtkCellArray>::New();
	int ist = 0;
	for (auto c: cont){
		vtkSmartPointer<vtkPolygon> aPolygon = vtkSmartPointer<vtkPolygon>::New();
		for (int i=0; i<c.n_points(); ++i) aPolygon->GetPointIds()->InsertNextId(ist++);
		vtkBndCell->InsertNextCell(aPolygon);
	}
	
	//4) Create a polydata to store the boundary.
	vtkSmartPointer<vtkPolyData> boundary = vtkSmartPointer<vtkPolyData>::New();
	boundary->SetPoints(vtkPoly->GetPoints());
	boundary->SetPolys(vtkBndCell);
	
	//5) perform Delaunay
	vtkSmartPointer<vtkDelaunay2D> delaunay = vtkSmartPointer<vtkDelaunay2D>::New();
	delaunay->SetInputData(vtkPoly);
	delaunay->SetSourceData(boundary);
	delaunay->Update();

	//6) build return vector from indicies of resulting cells
	vector<int> ret;
	auto pd = delaunay->GetOutput();
	for (int i=0; i<pd->GetNumberOfCells(); ++i){
		auto c = pd->GetCell(i);
		for (int j=0; j<3; ++j){
			ret.push_back(c->GetPointId(j));
		}
	}

	return ret;
}

void TriGrid::build(const vector<PContour>& cont, const shp_vector<Point>& pts){
	//triangulate
	auto tri = perform_triangulation(cont, pts);
	//build points
	for (auto c: cont){
		for (int i=0; i<c.n_points(); ++i){
			aa::add_shared(points, GPoint(*c.get_point(i)));
		}
	}
	for (auto p: pts) aa::add_shared(points, GPoint(*p));
	//build cells
	auto it = tri.begin();
	while (it!=tri.end()){
		aa::add_shared(cells, Cell());
		add_point_to_cell(cells.back().get(), points[*it++].get());
		add_point_to_cell(cells.back().get(), points[*it++].get());
		add_point_to_cell(cells.back().get(), points[*it++].get());
	}
	//indicies
	set_indicies();

	//force positive triangles (necessary)
	force_cells_ordering();

}

TriGrid::TriGrid(const vector<PContour>& cont, const shp_vector<Point>& pts): GridGeom(){
	build(cont, pts);
}

TriGrid::TriGrid(const vector<PContour>& cont, const vector<Point>& pts): GridGeom(){
	shp_vector<Point> p2;
	for(auto p: pts) aa::add_shared(p2, p);
	build(cont, p2);
}

std::set<Edge>& TriGrid::edges() const{
	if (_edges.size()==0) _edges = get_edges();
	return _edges;
}
std::map<GPoint*, vector<GPoint*>>& TriGrid::nodenodeI() const{
	if (_nodenodeI.size()==0){
		//1) get nodenode
		auto& eds = edges();
		for (auto& e: eds){
			auto it1 = _nodenodeI.emplace(points[e.p1].get(), vector<GPoint*>());
			it1.first->second.push_back(points[e.p2].get());
			auto it2 = _nodenodeI.emplace(points[e.p2].get(), vector<GPoint*>());
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

TriGrid TriGrid::refine_grid(double density){
	//1 get contours
	std::vector<PContour> bg_cont = get_contours();

	//2 define characteristic distance for each point in a contour vector
	vector<double> char_dist;
	for (auto& c2: bg_cont){
		vector<double> _chr = c2.chdist();
		std::copy(_chr.begin(), _chr.end(), std::back_inserter(char_dist));
	}

	//3 find refinement points on each grid edge
	shp_vector<Point> ref_pnt = ref_points(char_dist, density);
	
	//4 build another grid using refinement points
	TriGrid g3ref(bg_cont, ref_pnt);
	
	//5 smooth
	for (int i=0; i<3; ++i){ 
		for (int j=0; j<10; ++j) g3ref.smooth(1.0);
		vector<Point> p3;
		for (int j=0; j<g3ref.n_points(); ++j) p3.push_back(*g3ref.get_point(j));
		g3ref = TriGrid(bg_cont, p3);
	}

	//6 return
	return g3ref;
}

