#include "procgrid.h"

GridGeom GGeom::Constructor::RectGrid(const vector<double>& part_x,vector<double>& part_y){
	int Npts = part_x.size()*part_y.size();
	int Ncls = (part_x.size()-1)*(part_y.size()-1);
	vector<double> points;
	for (int j=0; j<part_y.size(); ++j){
		for (int i=0; i<part_x.size(); ++i){
			points.push_back(part_x[i]);
			points.push_back(part_y[j]);
		}
	}
	vector<int> cells;
	for (int j=0; j<part_y.size() - 1; ++j){
		for (int i=0; i<part_x.size() - 1; ++i){
			int i0 = j*part_x.size() + i;
			int i1 = j*part_x.size() + i + 1;
			int i2 = (j + 1) * part_x.size() + i + 1;
			int i3 = (j + 1) * part_x.size() + i;
			cells.push_back(4);
			cells.push_back(i0);
			cells.push_back(i1);
			cells.push_back(i2);
			cells.push_back(i3);
		}
	}
	return GridGeom(Npts, Ncls, &points[0], &cells[0]);
}

GridGeom GGeom::Constructor::RectGrid(Point p0, Point p1, int Nx, int Ny){
	vector<double> x(Nx), y(Ny);
	for (int i=0; i<x.size(); ++i){ x[i] = (p1.x-p0.x)*(double)i/(Nx-1)+p0.x; }
	for (int i=0; i<y.size(); ++i){ y[i] = (p1.y-p0.y)*(double)i/(Ny-1)+p0.y; }
	return RectGrid(x, y);
}

GridGeom GGeom::Constructor::RectGrid01(int Nx, int Ny){
	vector<double> x(Nx), y(Ny);
	for (int i=0; i<x.size(); ++i){ x[i] = (double)i/(Nx-1); }
	for (int i=0; i<y.size(); ++i){ y[i] = (double)i/(Ny-1); }
	return RectGrid(x, y);
}

GridGeom GGeom::Constructor::Ring(Point p0, double rad1, double rad2, int narc, int nrad){
	vector<double> pts;
	vector<int> cells;
	//1) build points
	for (int i=0; i<nrad+1; ++i){
		double r = rad1 + (rad2 - rad1)/nrad*i;
		for (int j=0; j<narc; ++j){
			double phi = 2*M_PI/narc*j;
			pts.push_back(r*cos(phi));
			pts.push_back(r*sin(phi));
		}
	}
	//2) build cells
	for (int i=0; i<nrad; ++i){
		for (int j=0; j<narc; ++j){
			int jnext = (j==narc-1)?0:j+1;
			cells.push_back(4);
			cells.push_back(narc*i+j);
			cells.push_back(narc*i+jnext);
			cells.push_back(narc*(i+1)+jnext);
			cells.push_back(narc*(i+1)+j);
		}
	}
	//3) build grid
	auto ret = GridGeom(pts.size()/2, cells.size()/5, &pts[0], &cells[0]);
	return ret;
}

GridGeom GGeom::Constructor::DeepCopy(const GridGeom& from){
	GridGeom ret;
	GGeom::Modify::DeepAdd(&from, &ret);
	return ret;
}

void GGeom::Modify::RemoveCells(GridGeom& grid, const std::vector<const Cell*>& cls){
	vector<int> bc(cls.size());
	std::transform(cls.begin(), cls.end(), bc.begin(), [](const Cell* c){ return c->get_ind(); });
	grid.remove_cells(bc);
}

void GGeom::Modify::AddCell(GridGeom& grid, const std::vector<Point>& cell){
	shared_ptr<Cell> c(new Cell);
	for (auto& p: cell) c->points.push_back(aa::add_shared(grid.points, GridPoint(p)));
	grid.cells.push_back(c);
	grid.set_indicies();
}

void GGeom::Modify::PointModify(GridGeom& grid, std::function<void(GridPoint*)> fun){
	std::for_each(grid.points.begin(), grid.points.end(),
			[&fun](shared_ptr<GridPoint> gp){ fun(gp.get()); });
}
void GGeom::Modify::CellModify(GridGeom& grid, std::function<void(Cell*)> fun){
	std::for_each(grid.cells.begin(), grid.cells.end(),
			[&fun](shared_ptr<Cell> gp){ fun(gp.get()); });
}

void GGeom::Modify::ShallowAdd(const GridGeom* from, GridGeom* to){
	std::copy(from->points.begin(), from->points.end(), std::back_inserter(to->points));
	std::copy(from->cells.begin(), from->cells.end(), std::back_inserter(to->cells));
	to->set_indicies();
}
void GGeom::Modify::ShallowAdd(const GridGeom* from, GridGeom* to, const std::vector<int>& icells, bool to_ind){
	std::set<int> usedpts;
	for (auto ic: icells){
		to->cells.push_back(from->cells[ic]);
		for (int ip=0; ip<from->cells[ic]->dim(); ++ip){
			usedpts.insert(from->cells[ic]->points[ip]->get_ind());
		}
	}
	for (auto ip:usedpts) to->points.push_back(from->points[ip]);
	if (to_ind) to->set_indicies();
}

void GGeom::Modify::ShallowAdd(const ShpVector<GridPoint>& from, GridGeom* to){
	std::copy(from.begin(), from.end(), std::back_inserter(to->points));
	to->set_indicies();
}

void GGeom::Modify::DeepAdd(const GridGeom* from, GridGeom* to){
	to->add_data(*from);
}

void GGeom::Repair::Heal(GridGeom& grid){
	grid.merge_congruent_points();
	grid.delete_unused_points();
	grid.force_cells_ordering();
}

bool GGeom::Repair::HasSelfIntersections(const GridGeom& grid){
	int n = 20;
	GGeom::Info::CellFinder cfinder(&grid, n, n);
	double ksieta[2];
	for (int i=0; i<n*n; ++i){
		auto cls = cfinder.CellsBySquare(i);
		//assemble cell edges
		std::set<Edge> edges;
		for (auto c: cls){
			int jprev = c->dim() - 1;
			for (int j=0; j<c->dim(); ++j){
				int i1 = c->points[jprev]->get_ind();
				int i2 = c->points[j]->get_ind();
				edges.emplace(i1, i2);
			}
		}
		//find intersections
		auto it1 = edges.begin();
		while (it1!= edges.end()){
			auto it2=std::next(it1);
			while (it2!=edges.end()){
				if (it1->p1 != it2->p1 &&
						it1->p1 != it2->p2 &&
						it1->p2 != it2->p1 &&
						it1->p2 != it2->p2){
					auto p1 = grid.get_point(it1->p1);
					auto p2 = grid.get_point(it1->p2);
					auto p3 = grid.get_point(it2->p1);
					auto p4 = grid.get_point(it2->p2);
					SectCross(*p1, *p2, *p3, *p4, ksieta);
					if (ksieta[0] > geps && ksieta[0] < 1-geps &&
						ksieta[1]>geps && ksieta[1] < 1-geps){
						return true;
					}
				}
				++it2;
			}
			++it1;
		}
	}
	return false;
}


ShpVector<GridPoint> GGeom::Info::SharePoints(const GridGeom& grid, const vector<int>& indicies){
	ShpVector<GridPoint> ret;
	for (int i: indicies) ret.push_back(grid.points[i]);
	return ret;
}

ShpVector<GridPoint> GGeom::Info::BoundaryPoints(const GridGeom& grid){
	auto pntset = grid.get_bnd_points();
	ShpVector<GridPoint> ret;
	for (auto pv: pntset) ret.push_back(grid.points[pv->get_ind()]);
	return ret;
}


HMCont2D::ContourTree GGeom::Info::Contour(const GridGeom& grid){
	//build subgrids.
	auto sg = Modify::SubGrids(grid);
	std::list<HMCont2D::Contour> all_contours;
	for (auto& g: sg){
		//indicies of 'g' cells/points are from 'grid'.
		//We need to use correct indicies for g in order to build edges.
		//At the end of this procedure we will bring all grid indicies back
		g.set_indicies();
		std::set<Edge> edges = g.get_edges();
		HMCont2D::ECollection ecol;
		for (auto& e: edges) if (e.is_boundary()){
			assert(g.points.size() > std::max(e.p1, e.p2));
			HMCont2D::Edge newe(g.points[e.p1].get(), g.points[e.p2].get());
			ecol.add_value(newe);
		}
		auto ac = HMCont2D::Assembler::AllContours(ecol);
		all_contours.insert(all_contours.end(), ac.begin(), ac.end());
	}
	all_contours.remove_if([](const HMCont2D::Contour& c){ return c.size()<3 || !c.is_closed(); });
	HMCont2D::ContourTree tree;
	for (auto& c: all_contours)  tree.AddContour(c); 
	//return correct indicies to grid
	const_cast<GridGeom*>(&grid)->set_indicies();
	return tree;
}

HMCont2D::Contour GGeom::Info::Contour1(const GridGeom& grid){
	HMCont2D::ContourTree ct = Contour(grid);
	return *ct.nodes[0];
}

namespace{
//recursive algorithm of connected cells traversal
void subgrids_add_cell(int i, vector<int>& v, const vector<vector<int>>& cell_cell, vector<int>& cind){
	cind[i]=1;
	v.push_back(i);
	for (auto ci: cell_cell[i]){
		if (cind[ci]==0) subgrids_add_cell(ci, v, cell_cell, cind);
	}
}
}//namespace
vector<GridGeom> GGeom::Modify::SubGrids(const GridGeom& grid){
	//cell->cell connectivity
	std::vector<std::vector<int>> cc = grid.cell_cell();
	//vector of cells use status (0 - unused, 1 - used)
	std::vector<int> cind(grid.n_cells(), 0);
	//build sets of single connected cells
	vector<vector<int>> sc_cells;
	while (1){
		auto fnd0 = std::find(cind.begin(), cind.end(), 0);
		//stop if all cells are already used
		if (fnd0 == cind.end()) break;
		sc_cells.push_back(vector<int>());
		subgrids_add_cell(fnd0-cind.begin(), sc_cells.back(), cc, cind);
	}

	//assemble new grids
	vector<GridGeom> ret;
	for (auto& sc: sc_cells){
		ret.push_back(GridGeom());
		auto& r = ret.back();
		ShallowAdd(&grid, &r, sc, false);
	}
	return ret;
}

void GGeom::Modify::SnapToContour(GridGeom& grid, const HMCont2D::Contour& cont,
		const std::vector<GridPoint*>& snap_nodes){
	//snapping nodes
	for (auto p: snap_nodes){
		//try to search amoung vertices
		Point* fpnt = HMCont2D::ECollection::FindClosestNode(cont, *p);
		if (Point::meas(*fpnt, *p)<geps*geps){
			p->set(*fpnt);
			continue;
		}
		//snap to edge
		auto fed = HMCont2D::ECollection::FindClosestEdge(cont, *p);
		HMCont2D::Edge* e = std::get<0>(fed);
		double w = std::get<2>(fed);
		p->set(Point::Weigh(*e->pstart, *e->pend, w));
	}

	auto gridbnd = GGeom::Info::Contour(grid);
	//all contour significant points weights
	std::map<double, Point*> contw;
	auto cp = cont.corner_points();
	for (auto p: cp){
		auto coord = cont.coord_at(*p);
		contw[std::get<1>(coord)] = p;
	}
	//copy one more time with w+1 for closed contours
	if (cont.is_closed()){
		auto it = contw.rbegin();
		while (it!=contw.rend()){
			contw[it->first + 1.0] = it->second;
			++it;
		}
	}

	//grid boundary points which lie on cont weights
	std::map<const GridPoint*, double> bndw;
	for (auto p: grid.get_bnd_points()){
		auto coord = cont.coord_at(*p);
		if (std::get<4>(coord)<geps) bndw[p] = std::get<1>(coord);
	}

	//assemble set of boundary grid edges
	std::set<Edge> bndeds;
	for (auto& e: grid.get_edges()) if (e.is_boundary()) bndeds.insert(e);
	//get all cont nodes between p1 and p2
	std::function<std::list<Point*>(double, double)> get_pts_between;
	get_pts_between = [&](double w1, double w2)->std::list<Point*>{
			std::list<Point*> ret;
			if (ISZERO(w2-w1)) return ret;
			if (w2<w1) w2+=1.0;
			assert(w1 < w2 && w2<2.0);
			for (auto it=contw.lower_bound(w1); it!=contw.end(); ++it){
				if (it->first>w1+geps && it->first<w2-geps){
					ret.push_back(it->second);
				}
				if (it->first>w2-geps) break;
			}
			return ret;
		};
	for (auto e: bndeds){
		const GridPoint *p1 = grid.get_point(e.p1), *p2 = grid.get_point(e.p2);
		if (e.cell_left < 0) std::swap(p1, p2);
		auto fnd1 = bndw.find(p1), fnd2 = bndw.find(p2);
		//if edge is not on contour ignore it
		if (fnd1==bndw.end() || fnd2==bndw.end()) continue;
		double w1 = fnd1->second, w2 = fnd2->second;
		std::list<Point*> ap = get_pts_between(w1, w2);
		//if no points between edge points ignore edge
		if (ap.size() == 0) continue;
		//add points to grid cell
		Cell* cell = grid.cells[e.any_cell()].get();
		auto i1 = std::find(cell->points.begin(), cell->points.end(), p1);
		auto i2 = std::find(cell->points.begin(), cell->points.end(), p2);
		assert(i2-i1==1 || (i2==cell->points.begin() && i1==cell->points.end()-1)); 
		//insert additional edges
		vector<GridPoint*> padds;
		for (auto it=ap.begin(); it!=ap.end(); ++it){
			grid.points.push_back(std::make_shared<GridPoint>(**it));
			padds.push_back(grid.points.back().get());
		}
		cell->points.insert(i1+1, padds.begin(), padds.end());
	}
	grid.set_indicies();
}

void GGeom::Modify::SimplifyBoundary(GridGeom& grid, double angle){
	assert(angle>-geps && angle<M_PI+geps);
	//get rid of nodes which are hanging + non-significant + boundary
	vector<std::map<int, int>> used_points_of_cells(grid.n_cells());
	std::map<int, std::set<int>> hangnodes;
	//function which define whether p2 in (p1, p2, p3) list is worth of deletion
	std::function<bool(const Point&,const Point&,const Point&)> isok;
	if (angle>M_PI-geps) isok = [](const Point&, const Point&, const Point&){return true;};
	else if (angle<geps) isok = [](const Point& p1, const Point& p2, const Point& p3)->bool{
		return fabs(triarea(p1, p2, p3)) < geps*geps;};
	else isok = [&](const Point& p1, const Point& p2, const Point& p3)->bool{
			double a = Angle(p1, p2, p3);
			if (a>M_PI) a=2*M_PI-a;
			a = fabs(M_PI-a);
			return a<angle;
		};
	auto addhn = [&](int ci, int gp, int ga){
		Cell* c = grid.cells[ci].get();
		hangnodes.emplace(ci, std::set<int>());
		int locp = 0, loca = 0, maxl = c->dim();
		while (c->points[locp]->get_ind()!=gp && locp<maxl) ++locp;
		while (c->points[loca]->get_ind()!=ga && loca<maxl) ++loca;
		assert(locp != maxl && loca != maxl);
		auto r = used_points_of_cells[ci].emplace(locp, loca);
		if (!r.second){
			int locb = r.first->second;
			if (isok(*c->points[loca],*c->points[locp], *c->points[locb]))
				hangnodes[ci].insert(locp);
		}
	};
	for (auto& e: grid.get_edges()) if (e.is_boundary()){
		int ei = e.any_cell();
		addhn(ei, e.p1, e.p2);
		addhn(ei, e.p2, e.p1);
	}
	auto it=hangnodes.begin();
	while (it!=hangnodes.end()){
		if (it->second.size() == 0) it = hangnodes.erase(it);
		else ++it;
	}
	if (hangnodes.size() > 0){
		for (auto& hn: hangnodes){
			//only if resulting cell will not become degenerate
			Cell cpcell(*grid.cells[hn.first]);
			aa::remove_entries(cpcell.points, hn.second);
			if (cpcell.dim()>2 && !cpcell.has_self_crosses())
				grid.cells[hn.first]->points = cpcell.points;
		}
		grid.delete_unused_points();
		grid.set_indicies();
	}
}

HMCont2D::Contour GGeom::Info::CellContour(const GridGeom& grid, int cell_ind){
	const Cell* c=grid.cells[cell_ind].get();
	vector<Point*> p(c->points.begin(), c->points.end());
	return HMCont2D::Constructor::ContourFromPoints(p, true);
}

BoundingBox GGeom::Info::BBox(const GridGeom& grid, double eps){
	auto ret = BoundingBox::Build(grid.points.begin(), grid.points.end());
	ret.widen(eps);
	return ret;
}

vector<double> GGeom::Info::Skewness(const GridGeom& grid){
	vector<double> ret(grid.n_cells());
	for (int i=0; i<grid.n_cells(); ++i){
		const Cell& c = *grid.get_cell(i);
		if (c.dim() < 3){
			ret[i] = 1.0;  //very bad cell anyway
			continue;
		}
		vector<double> angles(c.dim());
		for (int j=1; j<c.dim(); ++j){
			const Point& p0 = *c.get_point(j-1);
			const Point& p1 = *c.get_point(j);
			const Point& p2 = *c.get_point(j+1);
			angles[j] = Angle(p0, p1, p2);
		}
		angles[0] = M_PI *(c.dim() - 2) - 
			std::accumulate(angles.begin() + 1, angles.begin() + c.dim(), 0.0);
		auto minmax = std::minmax_element(angles.begin(), angles.end());
		double minv = *minmax.first;
		double maxv = *minmax.second;
		double refan = M_PI * (c.dim() - 2) / c.dim();
		ret[i] = std::max( (maxv-refan)/(M_PI-refan), (refan-minv)/refan );
		if (ret[i] > 1.0) ret[i] = 1.0;   //for non-convex cells
	}
	return ret;
}

double GGeom::Info::Area(const GridGeom& grid){
	return grid.area();
}

vector<double> GGeom::Info::CellAreas(const GridGeom& grid){
	vector<double> ret(grid.n_cells());
	for (int i=0; i<grid.n_cells(); ++i){
		ret[i] = grid.get_cell(i)->area();
	}
	return ret;
}

std::set<int> GGeom::Info::CellFinder::IndSet(const BoundingBox& bbox) const{
	Point b0 = bbox.BottomLeft() - p0, b1 = bbox.TopRight() - p0;
	int ix0 = std::floor(b0.x / Hx);
	int iy0 = std::floor(b0.y / Hy);
	int ix1 = std::floor(b1.x / Hx);
	int iy1 = std::floor(b1.y / Hy);

	std::set<int> ret;
	for (int j=iy0; j<iy1+1; ++j){
		if (j>=Ny) break;
		for (int i=ix0; i<ix1+1; ++i){
			if (i>=Nx) break;
			ret.insert(i+j*Nx);
		}
	}
	return ret;
}
GGeom::Info::CellFinder::CellFinder(const GridGeom* g, int nx, int ny):
		grid(g), Nx(nx), Ny(ny){
	BoundingBox bbox = GGeom::Info::BBox(*g);
	p0 = bbox.BottomLeft();
	Hx = bbox.lenx() / nx;
	Hy = bbox.leny() / nx;
	cells_by_square.resize(Nx*Ny);
	//associate cells with squares
	for (int i=0; i<g->n_cells(); ++i) AddCell(g->get_cell(i));
}

int GGeom::Info::CellFinder::GetSquare(const Point& p) const{
	int ix = std::floor((p.x - p0.x) / Hx);
	int iy = std::floor((p.y - p0.y) / Hy);
	if (ix >= Nx) ix = Nx-1; if (ix<0) ix = 0;
	if (iy >= Ny) iy = Ny-1; if (iy<0) iy = 0;
	return ix + iy*Nx;
}

const Cell* GGeom::Info::CellFinder::Find(const Point& p) const{
	auto candidates = CellCandidates(p);
	bool hint = true;
	for (auto c: candidates){
		int w = c->get_contour().is_inside(p, &hint);
		if (w != OUTSIDE) return c;
	}
	throw GGeom::EOutOfArea(p);
}

const Cell* GGeom::Info::CellFinder::FindExcept(const Point& p,
		const std::set<const Cell*>& exc) const{
	auto candidates = CellCandidates(p);
	bool hint = true;
	bool foundany = false;
	for (auto c: candidates){
		int w = c->get_contour().is_inside(p, &hint);
		if (w != OUTSIDE){
			foundany = true;
			if (exc.find(c) == exc.end()) return c;
		}
	}

	if (!foundany) throw GGeom::EOutOfArea(p);
	else return 0;
}

void GGeom::Info::CellFinder::AddCell(const Cell* c){
	auto bbox = BoundingBox::Build(c->points.begin(), c->points.end());
	bbox.widen(geps);
	std::set<int> inds = IndSet(bbox);
	for (auto j: inds){
		cells_by_square[j].push_back(c);
	}
}

const vector<const Cell*>& GGeom::Info::CellFinder::CellCandidates(const Point& p) const{
	//get square of point
	int ind = GetSquare(p);
	//return
	return cells_by_square[ind];
}

