#include "procgrid.h"
#include <stack>

GridGeom GGeom::Constructor::EmptyGrid(){
	return GridGeom();
}

GridGeom GGeom::Constructor::RectGrid(const vector<double>& part_x, const vector<double>& part_y){
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
	vector<double> x(Nx+1), y(Ny+1);
	for (int i=0; i<x.size(); ++i){ x[i] = (p1.x-p0.x)*(double)i/(Nx)+p0.x; }
	for (int i=0; i<y.size(); ++i){ y[i] = (p1.y-p0.y)*(double)i/(Ny)+p0.y; }
	return RectGrid(x, y);
}

GridGeom GGeom::Constructor::RectGrid01(int Nx, int Ny){
	vector<double> x(Nx+1), y(Ny+1);
	for (int i=0; i<x.size(); ++i){ x[i] = (double)i/(Nx); }
	for (int i=0; i<y.size(); ++i){ y[i] = (double)i/(Ny); }
	return RectGrid(x, y);
}

GridGeom GGeom::Constructor::Ring(Point p0, double rad1, double rad2, int narc, int nrad){
	assert(rad1>rad2);
	vector<double> pts;
	vector<int> cells;
	//1) build points
	for (int i=0; i<nrad+1; ++i){
		double r = rad1 + (rad2 - rad1)/nrad*i;
		for (int j=0; j<narc; ++j){
			double phi = 2*M_PI/narc*j;
			pts.push_back(r*cos(phi) + p0.x);
			pts.push_back(r*sin(phi) + p0.y);
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

GridGeom GGeom::Constructor::Circle(Point p0, double rad, int narc, int nrad, bool tri_center){
	double rad2 = rad*(1.0)/nrad;
	GridGeom ret = GGeom::Constructor::Ring(p0, rad, rad2, narc, nrad-1);
	if (!tri_center){
		auto newcell = aa::add_shared(ret.cells, Cell());
		for (int i=0; i<narc; ++i){
			int ind = ret.n_points() - narc + i;
			ret.add_point_to_cell(newcell, ret.get_point(ind));
		}
	} else {
		aa::add_shared(ret.points, GridPoint(p0));
		for (int i=0; i<narc; ++i){
			int ind1 = ret.n_points() - narc - 1 + i;
			int ind2 = (i == narc - 1) ? ret.n_points() - narc - 1
			                           : ind1 + 1;
			auto newcell = aa::add_shared(ret.cells, Cell());
			ret.add_point_to_cell(newcell, ret.get_point(ind1));
			ret.add_point_to_cell(newcell, ret.get_point(ind2));
			ret.add_point_to_cell(newcell, ret.get_point(ret.n_points() - 1));
		}
	}
	ret.set_indicies();
	return ret;
}

GridGeom GGeom::Constructor::DeepCopy(const GridGeom& from){
	GridGeom ret;
	GGeom::Modify::DeepAdd(&from, &ret);
	return ret;
}

GridGeom GGeom::Constructor::FromData(const ShpVector<GridPoint>& points, const ShpVector<Cell>& cells){
	GridGeom ret;
	ret.points = points;
	ret.cells = cells;
	ret.set_indicies();
	return ret;
}

GridGeom GGeom::Constructor::ExtractCells(const GridGeom& g, const std::vector<int>& cind, int policy){
	GridGeom ret;
	for (int i: cind) ret.cells.push_back(g.cells[i]);
	for (auto c: ret.cells){
		for (auto p: c->points){
			int ind = p->get_ind();
			if (ind >=0 && ind<g.n_points() && g.points[ind].get() == p){
				ret.points.push_back(g.points[ind]);
			} else {
				auto fnd = std::find_if(g.points.begin(), g.points.end(),
					[&p](shared_ptr<Point> sp){ return sp.get() == p; });
				if (fnd == g.points.end()) throw EException("cannot find point while exctracting cells");
				else ret.points.push_back(*fnd);
			}
		}
	}
	ret.points = aa::no_dublicates(ret.points);
	if (policy == 2) return ret;
	if (policy == 1){
		ret.set_indicies();
		return ret;
	}
	//backup old p indices
	std::vector<int> bu; bu.reserve(g.n_points());
	for (auto p: g.points) bu.push_back(p->get_ind());
	ret.set_indicies();
	ret = DeepCopy(ret);
	//set indicies back
	auto it = bu.begin();
	for (auto p: g.points) p->ind = *it++;
	//return deepcopied object
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

void GGeom::Modify::ClearAll(GridGeom& grid){
	grid.cells.clear();
	grid.points.clear();
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

void GGeom::Modify::ReallocatePrimitives(GridGeom& grid){
	ShpVector<GridPoint> newpoints(grid.n_points());
	ShpVector<Cell> newcells(grid.n_cells());

	for (int i=0; i<grid.n_points(); ++i){
		newpoints[i].reset(new GridPoint(*grid.get_point(i)));
	}
	for (int i=0; i<grid.n_cells(); ++i){
		newcells[i].reset(new Cell(*grid.get_cell(i)));
	}
	auto _indexer = aa::ptr_container_indexer(grid.points);
	_indexer.convert();
	for (int i=0; i<grid.n_cells(); ++i){
		for (int j=0; j<grid.cells[i]->dim(); ++j){
			int index = _indexer.index(grid.cells[i]->points[j]);
			newcells[i]->points[j] = newpoints[index].get();
		}
	}
	_indexer.restore();
	std::swap(grid.points, newpoints);
	std::swap(grid.cells, newcells);
	grid.set_indicies();
}

void GGeom::Repair::Heal(GridGeom& grid){
	grid.merge_congruent_points();
	grid.delete_unused_points();
	grid.force_cells_ordering();
}

void GGeom::Repair::RemoveShortEdges(GridGeom& grid, double ref_len){
	//1) calculate charachteristic cell sizes
	vector<double> csz(grid.n_cells());
	for (int i=0; i<grid.n_cells(); ++i){
		auto bbox = grid.cells[i]->bbox();
		csz[i] = ref_len*bbox.lendiag();
	}
	auto edges = grid.get_edges();
	//2) boundary points indicies
	std::set<int> bnd_points;
	for (auto e: edges) if (e.is_boundary()){
		bnd_points.insert(e.p1);
		bnd_points.insert(e.p2);
	}
	//3) collapse short edges
	for (auto e: edges){
		GridPoint& p1 = *grid.get_point(e.p1);
		GridPoint& p2 = *grid.get_point(e.p2);
		double elen = Point::dist(p1, p2);
		bool del = true;
		if (e.cell_left >=0 && elen > csz[e.cell_left]) del = false;
		if (e.cell_right >=0 && elen > csz[e.cell_right]) del = false;
		if (del){
			bool isbnd1 = bnd_points.find(p1.get_ind()) != bnd_points.end();
			bool isbnd2 = bnd_points.find(p2.get_ind()) != bnd_points.end();
			if (!isbnd1 && !isbnd2){
				double x = (p1.x + p2.x)/2.0, y = (p1.y + p2.y)/2.0;
				p1.set(x, y); p2.set(x, y);
			}
			else if (isbnd1 && !isbnd2) p2.set(p1);
			else if (!isbnd1 && isbnd2) p1.set(p2);
		}
	}
	
	//4) merge points
	Heal(grid);
}

bool GGeom::Repair::HasSelfIntersections(const GridGeom& grid){
	int n = 100;
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

void GGeom::Repair::NoConcaveCells(GridGeom& grid, double an){
	if (ISZERO(an)) an = 0;
	auto find_convex = [an](Cell& c)->int{
		for (int i=0; i<c.dim(); ++i){
			auto pm=c.get_point(i-1), p=c.get_point(i), pp=c.get_point(i+1);
			if (triarea(*pm, *p, *pp)>0) continue;
			double angle = Angle(*pm, *p, *pp) - M_PI;
			if (angle>=an) return i;
		}
		return -1;
	};
	auto has_cross = [](Cell& c, int i1, int i2)->bool{
		auto* p1 = c.get_point(i1);
		auto* p2 = c.get_point(i2);
		double ksi[2];
		for (int i=0; i<c.dim(); ++i){
			auto* p3 = c.get_point(i);
			auto* p4 = c.get_point(i+1);
			if (p3==p1 || p3==p2 || p4==p1 || p4==p2) continue;
			if (SectCross(*p1, *p2, *p3, *p4, ksi)) return true;
		}
		return false;
	};
	auto find_second_vert = [&](Cell& c, int v1)->int{
		int best = -1;
		double nrm=1e6;
		double ran = Angle(*c.get_point(v1-1), *c.get_point(v1), *c.get_point(v1+1))/2.0;
		for (int i=v1+2; i<=c.dim()+v1-2; ++i){
			if (has_cross(c, v1, i)) continue; //if edge doesn't cross cell
			double ian = Angle(*c.get_point(v1-1), *c.get_point(v1), *c.get_point(i));
			if (ian>=2*ran) continue;  //if edge is within the cell
			double diff = fabs(ran-ian);
			if (diff<nrm){nrm=diff; best=i;}
		}
		if (best<0) throw std::runtime_error("Failed to divide concave cell");
		return best;
	};
	auto divide_cell = [&](Cell& c, ShpVector<Cell>& add)->bool{
		int conv_vert = find_convex(c);
		if (conv_vert == -1) return false;
		int sec_vert = find_second_vert(c, conv_vert);
		if (conv_vert>sec_vert) std::swap(conv_vert, sec_vert);
		auto c1 = aa::add_shared(add, Cell());
		auto c2 = aa::add_shared(add, Cell());
		for (int i=conv_vert; i<=sec_vert; ++i){
			c1->points.push_back(const_cast<GridPoint*>(c.get_point(i)));
		}
		for (int i=sec_vert; i<=conv_vert+c.dim(); ++i){
			c2->points.push_back(const_cast<GridPoint*>(c.get_point(i)));
		}
		return true;
	};

	auto treat_cell = [&](Cell& c, ShpVector<Cell>& newc)->bool{
		ShpVector<Cell> arch;
		std::stack<Cell*> to_analyze;
		to_analyze.push(&c);
		while (to_analyze.size()>0){
			Cell* cur = to_analyze.top(); to_analyze.pop();
			bool divided=divide_cell(*cur, arch);
			if (divided){
				to_analyze.push(arch.end()[-1].get());
				to_analyze.push(arch.end()[-2].get());
			} else if (cur != &c) aa::add_shared(newc, Cell(*cur));
			else return false;
		}
		return true;
	};

	std::set<int> cells_to_delete;
	ShpVector<Cell> newcells;
	for (int i=0; i<grid.n_cells(); ++i){
		if (treat_cell(*grid.cells[i], newcells)) cells_to_delete.insert(i);
	}
	if (newcells.size() > 0){
		aa::remove_entries(grid.cells, cells_to_delete);
		std::copy(newcells.begin(), newcells.end(), std::back_inserter(grid.cells));
		grid.set_indicies();
	}
}

ShpVector<GridPoint> GGeom::Info::SharePoints(const GridGeom& grid){ return grid.points; }
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
	/*
	cind[i]=1;
	v.push_back(i);
	for (auto ci: cell_cell[i]){
		if (cind[ci]==0) subgrids_add_cell(ci, v, cell_cell, cind);
	}
	*/
	std::stack<int> cand; 
	cand.push(i);
	while (cand.size()>0){
		int j = cand.top(); cand.pop();
		if (cind[j] == 1) continue;
		cind[j] = 1;
		v.push_back(j);
		for (auto ci: cell_cell[j]){
			if (cind[ci] == 0) cand.push(ci);
		}
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
	if (sc_cells.size() == 1){
		ret.push_back(GridGeom());
		auto& r = ret.back();
		ShallowAdd(&grid, &r);
	} else for (auto& sc: sc_cells){
		ret.push_back(GridGeom());
		auto& r = ret.back();
		ShallowAdd(&grid, &r, sc, false);
	}
	return ret;
}

struct GGeom::Modify::_ShiftSnapPreCalc{
	GridGeom* g;
	HMCont2D::ContourTree gridbnd;
	std::map<double, Point*> contw;
	vector<Point*> cp;
	std::map<const GridPoint*, double> bndw;
	std::set<Edge> bndeds;

	_ShiftSnapPreCalc(GridGeom& grid, const HMCont2D::Contour& cont,
			const std::vector<GridPoint*>& snap_nodes): g(&grid){
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
		gridbnd = GGeom::Info::Contour(grid);
		//all contour significant points weights
		cp = cont.corner_points();
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
		for (auto p: grid.get_bnd_points()){
			auto coord = cont.coord_at(*p);
			if (std::get<4>(coord)<geps) bndw[p] = std::get<1>(coord);
		}
		//assemble set of boundary grid edges
		for (auto& e: grid.get_edges()) if (e.is_boundary()) bndeds.insert(e);
	}
	std::list<Point*> get_pts_between(double w1, double w2){
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
	}

	std::vector<std::tuple<GridPoint*, GridPoint*, std::list<Point*>, Cell*>>
	point_pairs(){
		std::vector<std::tuple<GridPoint*, GridPoint*, std::list<Point*>, Cell*>> ret;
		for (auto e: bndeds){
			GridPoint *p1 = g->points[e.p1].get(), *p2 = g->points[e.p2].get();
			if (e.cell_left < 0) std::swap(p1, p2);
			auto fnd1 = bndw.find(p1), fnd2 = bndw.find(p2);
			//if edge is not on contour ignore it
			if (fnd1==bndw.end() || fnd2==bndw.end()) continue;
			double w1 = fnd1->second, w2 = fnd2->second;
			std::list<Point*> ap = get_pts_between(w1, w2);
			if (ap.size() == 0) continue;
			Cell* cell = g->cells[e.any_cell()].get();
			ret.push_back(std::make_tuple(p1, p2, ap, cell));
		}
		return ret;
	}
};

void GGeom::Modify::SnapToContour(GridGeom& grid, const HMCont2D::Contour& cont,
		const std::vector<GridPoint*>& snap_nodes){
	auto proc = _ShiftSnapPreCalc(grid, cont, snap_nodes);
	for (auto& p: proc.point_pairs()){
		GridPoint* p1 = std::get<0>(p), *p2 = std::get<1>(p);
		auto& ap = std::get<2>(p);
		Cell* cell = std::get<3>(p);
		Cell c2(*cell);
		auto i1 = std::find(c2.points.begin(), c2.points.end(), p1);
		auto i2 = std::find(c2.points.begin(), c2.points.end(), p2);
		assert(i2-i1==1 || (i2==c2.points.begin() && i1==c2.points.end()-1)); 
		//try to insert additional edges
		ShpVector<GridPoint> padds;
		vector<GridPoint*> padds2;
		for (auto it=ap.begin(); it!=ap.end(); ++it){
			padds2.push_back(aa::add_shared(padds, GridPoint(**it)));
		}
		c2.points.insert(i1+1, padds2.begin(), padds2.end());
		if (!c2.has_self_crosses()){
			cell->points = c2.points;
			std::copy(padds.begin(), padds.end(), std::back_inserter(grid.points));
		}
	}
	grid.set_indicies();
}

namespace{
struct ShiftS{
	GridPoint* from;
	Point* to;
	double dist2;
};
bool operator<(const ShiftS& a, const ShiftS& b){ return a.dist2<b.dist2;}
}

void GGeom::Modify::ShiftToContour(GridGeom& grid, const HMCont2D::Contour& cont,
		const std::vector<GridPoint*>& snap_nodes){
	auto proc = _ShiftSnapPreCalc(grid, cont, snap_nodes);
	//get all possible shifts
	std::list<ShiftS> allshifts;
	for (auto p: proc.point_pairs()){
		GridPoint* p1 = std::get<0>(p), *p2 = std::get<1>(p);
		auto& ap = std::get<2>(p);
		ShiftS s1 {p1, *ap.begin(), Point::meas(*p1, **ap.begin())};
		ShiftS s2 {p2, *ap.rbegin(), Point::meas(*p2, **ap.rbegin())};
		allshifts.push_back(s1);
		allshifts.push_back(s2);
	}
	//leave only point on outer grid contour which are non-significant for cont
	for (auto& s: allshifts){
		for (auto p: proc.cp){
			if (Point::meas(*s.from, *p) < geps*geps){
				s.from = NULL;
				break;
			}
		}
	}
	allshifts.remove_if([](ShiftS s){ return s.from==NULL; });
	//sort allshifts
	allshifts.sort();
	//make shifts
	std::vector<std::vector<int>> point_cell = grid.point_cell();
	for (auto& s: allshifts) if (s.from != NULL){
		//try to shift
		Point bu(*s.from);
		s.from->set(*s.to);
		//check if all cells are still ok
		auto& ps = point_cell[s.from->get_ind()];
		bool ok = true;
		for (auto i=0; i<ps.size(); ++i){
			if (grid.cells[ps[i]]->has_self_crosses()){
				ok = false;
				break;
			}
		}
		if (ok){
			for (auto& s2: allshifts) if (&s2 != &s)
				if (s.to == s2.to || s.from == s2.from) s2.from = NULL;
		} else {
			s.from->set(bu);
		}
	}
	grid.set_indicies();
}

void GGeom::Modify::SnapAllBoundary(GridGeom& grid, const HMCont2D::ECollection& cont, int algo){
	auto gtree = GGeom::Info::Contour(grid);
	if (algo == 2){
		auto source = cont.all_points();
		for (auto p: gtree.all_points()){
			double minmeas = std::numeric_limits<double>::max();
			Point* to_p = 0;
			for (int i=0; i<source.size(); ++i){
				double m = Point::meas(*source[i], *p);
				if (m < minmeas){
					minmeas = m;
					to_p = source[i];
				}
				if (minmeas<geps*geps) break;
			}
			p->set(to_p->x, to_p->y);
		}
	} else {
		for (auto p: gtree.all_points()){
			auto closest = HMCont2D::ECollection::FindClosestEdge(cont, *p);
			if (std::get<1>(closest) > geps){
				HMCont2D::Edge* e = std::get<0>(closest);
				double w = std::get<2>(closest); 
				Point pnew = Point::Weigh(*e->pstart, *e->pend, w);
				p->set(pnew.x, pnew.y);
			}
		}
	}
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

bool GGeom::Info::Check(const GridGeom& grid){
	for (auto& c: grid.cells){
		double a = c->area();
		if (a<=0 || c->has_self_crosses()) return false;
	}
	return true;
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

