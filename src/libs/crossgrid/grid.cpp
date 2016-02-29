#include <assert.h>
#include <map>
#include "grid.h"
#include "addalgo.hpp"
#include "fileproc.h"
#include "trigrid.h"
#include "buffergrid.h"
#include "wireframegrid.h"
#include "procgrid.h"

void Edge::add_adj_cell(int cell_ind, int i1, int i2) const{
	if ((i1 == p1) && (i2 == p2)){
		cell_left = cell_ind;
	} else if ((i1==p2) && (i2==p1)){
		cell_right = cell_ind;
	} else {
		throw std::runtime_error("Invalid adjacent cell for an edge");
	}
}

int GridGeom::n_cellsdim() const{
	int ret = 0;
	for (auto c: cells){ ret+=c->dim(); }
	return ret;
}

PContour Cell::get_contour() const{
	return PContour::build(points);
}

double Cell::area() const{
	return PContour::build(points).area();
}

void Cell::check_ordering(){
	if (area()<0) std::reverse(points.begin(), points.end());
}

bool Cell::has_self_crosses() const{
	double ksi, ksieta[2];
	//1) points don't lie on edges
	for (int i=0; i<dim(); ++i){
		auto a = get_point(i);
		for (int j=1; j<dim()-1; ++j){
			auto b1 = get_point(i+j);
			auto b2 = get_point(i+j+1);
			isOnSection(*a, *b1, *b2, ksi);
			if (ksi>-geps && ksi<1+geps) return true;
		}
	}
	//2) edges don't cross each other
	for (int i=0; i<dim()-1; ++i){
		auto a1 = get_point(i);
		auto a2 = get_point(i+1);
		int lastp = (i==0) ? dim() - 1 : dim();
		for (int j=i+2; j<lastp; ++j){
			auto b1 = get_point(j);
			auto b2 = get_point(j+1);
			if (SectCross(*a1, *a2, *b1, *b2, ksieta)) return true;
		}
	}
	return false;
}

void GridGeom::add_data(const GridGeom& g){
	int start_p_index = n_points();
	for (auto p: g.points) aa::add_shared(points, GridPoint(*p));
	for (auto c: g.cells){
		auto nc = aa::add_shared(cells, Cell());
		for (int j=0; j<c->dim(); ++j){
			int pind = start_p_index + c->points[j]->ind;
			nc->points.push_back(points[pind].get());
		}
	}
	set_indicies();
}

void GridGeom::add_data(const GridGeom& g, const std::vector<int>& cls){
	bool need_merge = n_points()==0;
	if (cls.size() == g.n_cells()) add_data(g);
	else{
		std::map<int, GridPoint*> added_points;
		for (auto ci: cls){
			const Cell* c = g.get_cell(ci);
			Cell* newc = aa::add_shared(cells, Cell());
			for (int j=0; j<c->dim(); ++j){
				int orig_ind = c->get_point(j)->get_ind();
				auto fnd = added_points.find(orig_ind);
				GridPoint* cpnt = (fnd==added_points.end()) ? 
					aa::add_shared(points, GridPoint(*g.get_point(orig_ind))) :
					fnd->second;
				add_point_to_cell(newc, cpnt);
			}
		}
		set_indicies();
	}
	if (need_merge) merge_congruent_points();
}

void GridGeom::remove_cells(const vector<int>& bad_cells){
	if (bad_cells.size()==0) return;
	std::set<int> bc(bad_cells.begin(), bad_cells.end());
	aa::remove_entries(cells, bc);
	delete_unused_points();
	set_indicies();
}

void GridGeom::leave_only(const ContoursCollection& cnt){
	vector<Point> intp = cells_internal_points();
	vector<int> bad_cells = std::get<2>(cnt.filter_points_i(intp));
	remove_cells(bad_cells);
}
void GridGeom::exclude_area(const ContoursCollection& cnt){
	vector<Point> intp = cells_internal_points();
	vector<int> bad_cells = std::get<0>(cnt.filter_points_i(intp));
	remove_cells(bad_cells);
}

void GridGeom::merge_congruent_points(){
	if (n_points() == 0) return;
	std::map<int, int> cp;
	NodeFinder fnd(outer_rect());
	for (int i=0; i<n_points(); ++i){
		const GridPoint* p = get_point(i);
		auto pf = static_cast<const GridPoint*>(fnd.add(p));
		if (pf!=p) cp.emplace(i, pf->get_ind());
	}

	if (cp.size()>0){
		for (auto c: cells){
			for (int j=0; j<c->dim(); ++j){
				auto fnd = cp.find(c->get_point(j)->get_ind());
				if (fnd!=cp.end()){
					change_point_of_cell(c.get(), j, points[fnd->second].get());
				}
			}
		}
		delete_unused_points();
		set_indicies();
	}
}

GridGeom::GridGeom(const GridGeom& g){
	add_data(g);
}

GridGeom::GridGeom(GridGeom&& g){
	std::swap(points, g.points);
	std::swap(cells, g.cells);
}

void GridGeom::swap_data(GridGeom& g1, GridGeom& g2){
	std::swap(g1.points, g2.points);
	std::swap(g1.cells, g2.cells);
}

GridGeom& GridGeom::operator=(GridGeom g){
	std::swap(points, g.points);
	std::swap(cells, g.cells);
	return *this;
}

GridGeom::GridGeom(int Npts, int Ncells, double* pts, int* cls){
	//1) points
	for (int i=0; i<Npts; ++i){
		double x = *pts++;
		double y = *pts++;
		aa::add_shared(points, GridPoint(x,y,i));
	}
	//2) cells
	for (int i=0; i<Ncells; ++i){
		int n = *cls++;
		auto newcell = aa::add_shared(cells, Cell(i));
		for (int j=0; j<n; ++j){
			newcell->points.push_back(points[*cls++].get());
		}
	}
}
void GridGeom::set_indicies(){
	int i=0, j=0;
	for (auto p: points) p->ind = i++;
	for (auto c: cells)  c->ind = j++;
}

void GridGeom::delete_unused_points(){
	set_indicies();
	std::vector<int> usage(n_points(), 0);
	for (int i=0; i<n_cells(); ++i){
		for (int j=0; j<cells[i]->dim(); ++j){
			++usage[cells[i]->points[j]->ind];
		}
	}
	std::set<int> unused_points;
	for (size_t i=0; i<usage.size(); ++i){
		if (usage[i]==0){
			unused_points.insert(i);
		}
	}
	aa::remove_entries(points, unused_points);
}

ScaleBase GridGeom::do_scale(){
	return ScaleBase::p_doscale(points.begin(), points.end());
}
void GridGeom::do_scale(const ScaleBase& sc){
	sc.p_scale(points.begin(), points.end());
}
void GridGeom::undo_scale(const ScaleBase& sc){
	sc.p_unscale(points.begin(), points.end());
}

std::set<Edge> GridGeom::get_edges() const{
	std::set<Edge> ret;
	for (int i=0; i<n_cells(); ++i){
		int jprev = cells[i]->dim()-1;
		for (int j=0;j<cells[i]->dim();++j){
			int i1 = cells[i]->points[jprev]->ind;
			int i2 = cells[i]->points[j]->ind;
			auto ins = ret.emplace(i1, i2);
			ins.first->add_adj_cell(i, i1, i2);
			jprev = j;
		}
	}
	return ret;
}

std::set<const GridPoint*> GridGeom::get_bnd_points() const{
	std::set<const GridPoint*> ret;
	for (auto e: get_edges()) if (e.is_boundary()){
		ret.insert(points[e.p1].get());
		ret.insert(points[e.p2].get());
	}
	return ret;
}

std::vector<PContour> GridGeom::get_contours() const {
	//1) get boundary edges
	std::set<Edge> bnd_edges;
	for (auto e: get_edges()){
		if (e.is_boundary()){
			bnd_edges.insert(e);
		}
	}
	//2) connect edges into contours
	std::vector<PContour> ret;
	while (bnd_edges.size()>0){
		std::vector<int> cont;
		cont.push_back(bnd_edges.begin()->p1);
		cont.push_back(bnd_edges.begin()->p2);
		bool dir = bnd_edges.begin()->cell_right<0;
		int last_cell = bnd_edges.begin()->any_cell();
		bnd_edges.erase(bnd_edges.begin());
		//building a closed contour
		while (cont.front()!=cont.back()){
			std::set<Edge>::iterator e;
			std::vector<std::set<Edge>::iterator> candidates;
			for (e = bnd_edges.begin(); e!=bnd_edges.end(); ++e){
				if (e->p1 == cont.back() || e->p2 == cont.back())
					candidates.push_back(e);
			}
			//if more then 1 candidate prefer the one with the same element
			if (candidates.size() == 1){
				e = candidates[0];
			} else if (candidates.size() > 1){
				auto fnd = std::find_if(candidates.begin(), candidates.end(),
						[&last_cell](std::set<Edge>::iterator it){
							return it->any_cell() == last_cell;
						});
				if (fnd!=candidates.end()) e = *fnd;
				else e = candidates[0];
			} else { 
				throw std::runtime_error("Cannot detect a closed contour");
			}
			//add next point
			if (e->p1==cont.back()){
				cont.push_back(e->p2);
			} else cont.push_back(e->p1);
			last_cell = e->any_cell();
			//erase used edge
			bnd_edges.erase(e);
		}
		//remove last point which equals the first one
		cont.pop_back();
		//change a direction if point1->point2 has right cell
		if (!dir) std::reverse(cont.begin(), cont.end());
		//build a contour
		PContour c;
		for (auto i: cont) c.add_point(points[i].get());
		ret.push_back(c);
	}
	return ret;
}

std::pair<Point, Point> GridGeom::outer_rect() const{
	double x0(points[0]->x), y0(points[0]->y);
	double x1(points[0]->x), y1(points[0]->y);
	for (int i=1; i<n_points(); ++i){
		if (points[i]->x > x1) x1 = points[i]->x;
		if (points[i]->y > y1) y1 = points[i]->y;
		if (points[i]->x < x0) x0 = points[i]->x;
		if (points[i]->y < y0) y0 = points[i]->y;
	}
	return std::make_pair(Point(x0, y0), Point(x1, y1));
}


void GridGeom::force_cells_ordering(){
	for (auto c: cells) c->check_ordering();
}
std::vector<std::vector<int>> GridGeom::cell_cell() const{
	//cell->cell connectivity
	std::vector<std::vector<int>> cc(n_cells());
	auto edges = get_edges();
	for (auto& e: edges){
		if (e.cell_left>=0 && e.cell_right>=0){
			cc[e.cell_left].push_back(e.cell_right);
			cc[e.cell_right].push_back(e.cell_left);
		}
	}
	return cc;
}

namespace{
//recursive algorithm of connected cells traversal
void add_cell(int i, vector<int>& v, const vector<vector<int>>& cell_cell, vector<int>& cind){
	cind[i]=1;
	v.push_back(i);
	for (auto ci: cell_cell[i]){
		if (cind[ci]==0) add_cell(ci, v, cell_cell, cind);
	}
}
}//namespace

ShpVector<GridGeom> GridGeom::subdivide() const{
	auto cc = cell_cell();
	//vector of cells use status (0 - unused, 1 - used)
	std::vector<int> cind(n_cells(), 0);
	//build sets of single connected cells
	vector<vector<int>> sc_cells;
	while (1){
		auto fnd0 = std::find(cind.begin(), cind.end(), 0);
		//stop if all cells are already used
		if (fnd0 == cind.end()) break;
		sc_cells.push_back(vector<int>());
		add_cell(fnd0-cind.begin(), sc_cells.back(), cc, cind);
	}

	//assemble new grids
	ShpVector<GridGeom> ret;
	for (auto& sc: sc_cells){
		ret.push_back(std::shared_ptr<GridGeom>(new GridGeom()));
		auto r = ret.back().get();
		r->add_data(*this, sc);
	}
	return ret;
}

void GridGeom::no_excessive_bnodes(){
	//get rid of nodes which are hanging + non-significant + boundary
	vector<std::map<int, int>> used_points_of_cells(n_cells());
	std::map<int, std::set<int>> hangnodes;
	auto addhn = [&](int ci, int gp, int ga){
		Cell* c = cells[ci].get();
		hangnodes.emplace(ci, std::set<int>());
		int locp = 0, loca = 0, maxl = c->dim();
		while (c->points[locp]->get_ind()!=gp && locp<maxl) ++locp;
		while (c->points[loca]->get_ind()!=ga && loca<maxl) ++loca;
		assert(locp != maxl && loca != maxl);
		auto r = used_points_of_cells[ci].emplace(locp, loca);
		if (!r.second){
			int locb = r.first->second;
			if (ISZERO(triarea(*c->points[loca],
					*c->points[locb], *c->points[locp]))){
				hangnodes[ci].insert(locp);
			}
		}
	};
	for (auto& e: get_edges()) if (e.is_boundary()){
		int ei = e.any_cell();
		addhn(ei, e.p1, e.p2);
		addhn(ei, e.p2, e.p1);
	}
	if (hangnodes.size() > 0){
		for (auto& hn: hangnodes){
			aa::remove_entries(cells[hn.first]->points, hn.second);
		}
		delete_unused_points();
		set_indicies();
	}
}

GridGeom* GridGeom::combine(GridGeom* gmain, GridGeom* gsec){
	//1) build grids contours
	auto maincont = gmain->get_contours_collection();
	auto seccont = gsec->get_contours_collection();
	//2) input data to wireframe format
	PtsGraph wmain(*gmain);
	PtsGraph wsec(*gsec);
	//3) cut outer grid with inner grid contour
	wmain = PtsGraph::cut(wmain, seccont, INSIDE);
	//4) overlay grids
	wmain = PtsGraph::overlay(wmain, wsec);
	//5) build single connected grid
	GridGeom* ret = new GridGeom(wmain.togrid());
	//6) filter out all cells which lie outside gmain and gsec contours
	//if gmain and gsec are not simple structures
	if (maincont.n_cont()>1 || seccont.n_cont()>1){
		vector<Point> pts = ret->cells_internal_points();
		auto mainfilt = std::get<2>(maincont.filter_points_i(pts));
		auto secfilt = std::get<2>(seccont.filter_points_i(pts));
		vector<int> bad_cells;
		for (int i: mainfilt){
			if (std::find(secfilt.begin(), secfilt.end(), i)!=secfilt.end()){
				bad_cells.push_back(i);
			}
		}
		ret->remove_cells(bad_cells);
	}
	//7) get rid of hanging + non-significant + boundary nodes
	//   They appear if boundaries of gmain and gsec partly coincide
	ret->no_excessive_bnodes();
	return ret;
}

namespace { //cross_grids helper functions
vector<Point> intersection_points(const GridGeom& gmain, const GridGeom& gsec){
	//finds all points which lie on the intersection of gmain/gsec areas.
	//Only those points which lie in (gmain-gsec) area are considered
	//to exclude cases when intersection point is not actual after crossing
	vector<Point> ret;
	auto cmain = GGeom::Info::Contour(gmain);
	auto simp_main = HMCont2D::ContourTree::Simplified(cmain);
	auto csec = GGeom::Info::Contour(gsec);
	auto simp_sec = HMCont2D::ContourTree::Simplified(csec);
	auto cross = HMCont2D::Clip::Difference(simp_main, simp_sec);
	HMCont2D::Clip::Heal(cross);
	auto bbox1 = HMCont2D::ECollection::BBox(simp_sec, geps);
	for (auto p: cross.all_points()){
		if (bbox1.whereis(*p) == OUTSIDE) continue;
		auto mainfnd = HMCont2D::ECollection::FindClosestEdge(simp_main, *p);
		auto secfnd = HMCont2D::ECollection::FindClosestEdge(simp_sec, *p);
		if (ISZERO(std::get<1>(mainfnd)) && ISZERO(std::get<1>(secfnd))){
			ret.push_back(*p);
		}
	}
	return ret;
}
}

GridGeom* GridGeom::cross_grids(GridGeom* gmain_inp, GridGeom* gsec_inp, 
		double buffer_size, double density, bool preserve_bp, bool empty_holes,
		CrossGridCallback::func cb){
	//initial scaling before doing anything
	if (cb("Building grid cross", "Scaling", 0, -1) == CrossGridCallback::CANCEL) return 0;
	auto sc = gmain_inp->do_scale();
	gsec_inp->do_scale(sc);
	buffer_size/=sc.L;

	//procedure assigning
	GridGeom* gmain = gmain_inp;
	GridGeom* gsec = gsec_inp;
	std::shared_ptr<GridGeom> _gm, _gs;  //this will store modified gmain and gsec if needed
	auto cmain  = gmain->get_contours();
	auto csec  = gsec->get_contours();

	//1 ---- if secondary holes are empty -> remove secondary top level area from main grid
	if (cb("Building grid cross", "Boundary analyze", 0.05, -1) == CrossGridCallback::CANCEL) return 0;
	if (empty_holes){
		ContoursCollection c1(csec);
		if (c1.n_cont() != c1.n_inner_cont()){
			PointsContoursCollection c2(c1.cut_by_level(0, 0));
			_gm.reset(grid_minus_cont(gmain, &c2, true, cb));
			if (_gm.get() == 0) return 0;
			gmain = _gm.get();
			cmain = gmain->get_contours();
		}
	}
	
	if (cb("Building grid cross", "Boundary analyze", 0.20, -1) == CrossGridCallback::CANCEL) return 0;
	//2 ---- find contours intersection points and place gsec nodes there
	if (!preserve_bp && buffer_size>geps){
		//find all intersections
		vector<Point> bnd_intersections = intersection_points(*gmain, *gsec);
		/*
		for (int i=0; i<csec.size(); ++i){
			auto bint = csec[i].intersections(cmain);
			for (auto in: bint){
				if (fabs(in.first - round(in.first))>geps){
					bnd_intersections.push_back(in.second);
				}
			}
		}
		*/
		//move secondary grid boundary points to intersection points
		if (bnd_intersections.size() > 0){
			_gs.reset(new GridGeom(*gsec));
			gsec = _gs.get();
			gsec->move_boundary_points(bnd_intersections);
			csec  = gsec->get_contours();
		}
	}
	
	//3 ---- combine grids without using buffers
	if (cb("Building grid cross", "Combining", 0.25, -1) == CrossGridCallback::CANCEL) return 0;
	std::unique_ptr<GridGeom> comb(GridGeom::combine(gmain, gsec));

	//4 ---- fill buffer zone
	//zero buffer requires no further actions
	if (buffer_size>geps){
	
		//loop over each single connected grid
		if (cb("Building grid cross", "Subdivision", 0.45, -1) == CrossGridCallback::CANCEL) return 0;
		auto sg = comb->subdivide();

		//callback percentages 
		double cb_step = 1.0/(sg.size()*csec.size());
		double cb_cur = 0;
		for (auto grid: sg){
			//loop over each secondary contour
			for (auto c: csec){
				cb_cur+=cb_step;
				//1. filter out a grid from buffer zone for the contour
				if (cb("Building grid cross", "Filling buffer", 0.5, cb_cur-1.0*cb_step) 
						== CrossGridCallback::CANCEL) return 0;
				BufferGrid bg(*grid, c, buffer_size);

				//continue if buffergrid is empty
				if (bg.num_orig_cells()==0) continue; 

				//2. perform triangulation of buffer grid area
				if (cb("Building grid cross", "Filling buffer", 0.5, cb_cur-0.8*cb_step)
						== CrossGridCallback::CANCEL) return 0;
				auto bgcont = bg.boundary_info(preserve_bp);
				auto g3 = TriGrid::TriangulateArea(std::get<0>(bgcont), std::get<1>(bgcont), 1);

				//3. change the internal of bg by g3ref grid
				if (cb("Building grid cross", "Filling buffer", 0.5, cb_cur-0.6*cb_step)
						== CrossGridCallback::CANCEL) return 0;
				bg.change_internal(*g3);
				
				//4. update original grid using new filling of buffer grid
				if (cb("Building grid cross", "Filling buffer", 0.5, cb_cur-0.3*cb_step)
						== CrossGridCallback::CANCEL) return 0;
				bg.update_original();
			}
		}
		if (cb("Building grid cross", "Filling buffer", 0.5, 1.0)
						== CrossGridCallback::CANCEL) return 0;

		if (cb("Building grid cross", "Merging", 0.85, -1)
						== CrossGridCallback::CANCEL) return 0;
		//connect single connected grids
		comb->clear();
		for (size_t i=0; i<sg.size(); ++i){
			comb->add_data(*sg[i]);
		}
		if (sg.size()>0) comb->merge_congruent_points();
	}

	//5. --- scale back after all procedures have been done and return
	cb("Building grid cross", "Unscaling", 0.95, -1);
	gmain_inp->undo_scale(sc);
	gsec_inp->undo_scale(sc);
	comb->undo_scale(sc);
	comb->set_indicies();

	cb("Building grid cross", "Done", 1.0, -1);
	return comb.release();
}

void GridGeom::change_internal(const GridGeom& gg){
	//TODO leave pointers to boundary points
	points.clear();
	cells.clear();
	for (int i = 0; i<gg.n_points(); ++i){
		auto p = gg.get_point(i);
		aa::add_shared(points, GridPoint(p->x, p->y));
	}
	for (int i = 0; i<gg.n_cells(); ++i){
		auto newc = aa::add_shared(cells, Cell());
		auto oldc = gg.get_cell(i);
		for (auto j=0; j<oldc->dim(); ++j){
			int ind = oldc->get_point(j)->get_ind();
			add_point_to_cell(newc, points[ind].get());
		}
	}
	set_indicies();
}

vector<Point> GridGeom::cells_internal_points() const{
	vector<Point> ret; ret.reserve(n_cells());
	for (auto& c: cells){
		//1) find if convex and any positive triangle
		bool conv = true;
		int anypos = -1;
		for (int i=0; i<c->dim(); ++i){
			auto p1 = c->get_point(i-1);
			auto p2 = c->get_point(i);
			auto p3 = c->get_point(i+1);
			double ar = triarea(*p1, *p2, *p3);
			if (ar<-geps2){
				conv = false;
				break;
			} else if (ar>geps2){
				anypos = i;
			}
		}
		//2) if cell is convex: get center of any triangle
		if (conv && anypos>-1){
			auto p1 = c->get_point(anypos-1);
			auto p2 = c->get_point(anypos);
			auto p3 = c->get_point(anypos+1);
			ret.push_back((*p1 + *p2 + *p3)/3.0);
		} else if (!conv){
		//3) if cell is non-convex: ray crossing algo
			PContour cont;
			for (int i=0; i<c->dim();
				cont.add_point(const_cast<GridPoint*>(c->get_point(i++))));
			ret.push_back(cont.inside_point_ca());
		} else {
			throw std::runtime_error("singular cell");
		}
	}
	return ret;
}

GridGeom* GridGeom::grid_minus_cont(GridGeom* g, PointsContoursCollection* c,
		bool is_inner,
		CrossGridCallback::func cb){
	//initial scaling before doing anything
	const char* com = "Excluding contour from grid";
	if (cb(com, "Scaling", 0, -1) == CrossGridCallback::CANCEL) return 0;
	auto sc = g->do_scale();
	c->do_scale(sc);

	//1) build a grid boundary
	auto gbnd = g->get_contours_collection();

	//2) force is_inner=false by addition new contour to c
	ContoursCollection cp = c->shallow_copy();
	Contour bounding_cnt;
	if (is_inner==true){
		CGBoundingBox bbox2(cp);
		CGBoundingBox bbox1(gbnd);
		bounding_cnt = CGBoundingBox({bbox1, bbox2}, 1.0).get_contour();
		cp.add_contour(bounding_cnt);
	}

	//3) loop over each inner
	int cball = cp.n_inner_cont();
	int cbcur = 0;
	vector<GridGeom> gg;
	for (int i=0; i<cp.n_cont(); ++i) if (cp.is_inner(i)){
		double k = 0.05 + 0.90*(double)cbcur++/cball;
		if (cb(com, "Contours loop", k, 0) == CrossGridCallback::CANCEL) return 0;

		//leave only inner + outers in contourscollection
		ContoursCollection _ctmp = cp.level_01(i);

		//build graph
		PtsGraph graph(*g);

		// ------- add edges
		if (cb(com, "Contours loop", k, 0.1) == CrossGridCallback::CANCEL) return 0;
		//loop over each outer
		for (auto& child: _ctmp.get_childs(0)){
			graph.add_edges(*child);
		}
		if (cb(com, "Contours loop", k, 0.4) == CrossGridCallback::CANCEL) return 0;
		//treat inner by itself. always add contour (ignore_trivial=true)
		//because trivial edges will be the last to present in graph
		//after exclusion (not actual)
		graph.add_edges(*cp.get_contour(i));

		if (cb(com, "Contours loop", k, 0.5) == CrossGridCallback::CANCEL) return 0;
		// ------- exclude edges
		graph.exclude_area(gbnd,   OUTSIDE);
		graph.exclude_area(_ctmp,  OUTSIDE);

		gg.push_back(graph.togrid());

		if (cb(com, "Contours loop", k, 0.8) == CrossGridCallback::CANCEL) return 0;
		//remove elements which can present in the grid due
		//to outer contour exclusion
		if (_ctmp.num_childs(0) > 0) gg.back().leave_only(_ctmp);
		if (gbnd.n_cont() > 0) gg.back().leave_only(gbnd);
	}
	
	//4) assemble a grid from graph
	if (cb(com, "Assembling", 0.95, -1) == CrossGridCallback::CANCEL) return 0;
	auto ret = new GridGeom();
	for (auto& g: gg) ret->add_data(g);

	//5) unscaling
	g->undo_scale(sc);
	c->undo_scale(sc);
	ret->undo_scale(sc);

	cb(com, "Done", 1.0, -1);
	return ret;
}

vector<vector<int>> GridGeom::point_cell_tab() const{
	vector<vector<int>> ret(n_points());
	for (int i=0; i<n_cells(); ++i){
		auto c = get_cell(i);
		for (auto j=0; j<c->dim(); ++j){
			ret[c->get_point(j)->ind].push_back(c->ind);
		}
	}
	return ret;
}

void GridGeom::move_boundary_points(const vector<Point>& bp){
	auto bnd = get_contours();
	double ksi;
	auto pct = point_cell_tab();
	for (auto p: bp){
		//1) find edge in a loop for all boundary edges
		const GridPoint *b1=0, *b2=0;
		for (auto& c: bnd){
			for (int i=0; i<c.n_points(); ++i){
				auto p1 = c.get_point(i), p2 = c.get_point(i+1);
				if (isOnSection(p, *p1, *p2, ksi)){
					if (!c.is_corner_point(i))
						b1 = static_cast<const GridPoint*>(p1);
					if (!c.is_corner_point(i+1))
						b2 = static_cast<const GridPoint*>(p2);
					//no points which can be moved
					if (b1 == 0 && b2 == 0) goto edge_not_found;
					//sort: b1 should be better then b2
					if (b1!=0 && b2!=0 && ksi>0.5) std::swap(b1,b2);
					if (b1==0) std::swap(b1, b2);
					goto edge_found;
				}
			}
		}
		edge_not_found:
			continue;
		edge_found:
			//2) move point: try b1, if failed use b2
			if (!move_point(b1, p, pct) && b2!=0) move_point(b2, p, pct);
	}
}

bool GridGeom::move_point(const GridPoint* srcp, const Point& tarp, const vector<vector<int>>& point_cell){
	Point oldpoint(*srcp);
	GridPoint* s = points[srcp->ind].get();
	s->x = tarp.x; s->y = tarp.y;
	//check cells for self crosses and positivity
	bool ok=true;
	for (auto ci: point_cell[s->ind]){
		Cell* cell = cells[ci].get();
		if (cell->area() <= geps2 || cell->has_self_crosses()){
			ok = false;
			break;
		}
	}
	if (ok){
		return true;
	}else{
		s->x = oldpoint.x; s->y = oldpoint.y;
		return false;
	}
}

GridGeom GridGeom::sum(const std::vector<GridGeom>& g){
	if (g.size() == 0) return GridGeom();
	GridGeom ret(g[0]);
	for (int i=1; i<g.size(); ++i){
		ret.add_data(g[i]);
	}
	return ret;
}

GridGeom GridGeom::sum(const std::vector<GridGeom*>& g){
	if (g.size() == 0) return GridGeom();
	GridGeom ret(*g[0]);
	for (int i=1; i<g.size(); ++i){
		ret.add_data(*g[i]);
	}
	return ret;
}

double GridGeom::area() const{
	double a=0;
	for (int i=0; i<n_cells(); ++i) a+=get_cell(i)->area();
	return a;
}
