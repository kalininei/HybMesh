#include <assert.h>
#include <map>
#include "grid.h"
#include "addalgo.hpp"
#include "fileproc.h"
#include "trigrid.h"
#include "buffergrid.h"
#include "wireframegrid.h"
#include "bgeom.h"

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

double Cell::area() const{
	return PContour::build(points).area();
}
void Cell::check_ordering(){
	if (area()<0) std::reverse(points.begin(), points.end());
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

void GridGeom::merge_congruent_points(){
	std::map<int, int> cp;
	for (int i=0; i<n_points(); ++i){
		for (int j=i+1; j<n_points(); ++j){
			if (*get_point(i)==*get_point(j))
				cp.emplace(j,i);
		}
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
		bnd_edges.erase(bnd_edges.begin());
		//building a closed contour
		while (cont.front()!=cont.back()){
			std::set<Edge>::iterator e;
			for (e = bnd_edges.begin(); e!=bnd_edges.end(); ++e){
				if (e->p1==cont.back()){
					cont.push_back(e->p2);
					break;
				} else if (e->p2==cont.back()){
					cont.push_back(e->p1);
					break;
				}
			}
			if (e!=bnd_edges.end()){
				bnd_edges.erase(e);
			} else {
				throw std::runtime_error("Cannot detect a closed contour");
			}
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

GridGeom GridGeom::remove_area(const PContour& cont){
	//1) filter points which lies within the contour
	std::vector<Point*> p_all, p_inner, p_outer;
	for (auto& p: points) p_all.push_back(p.get());
	cont.select_points(p_all, p_inner, p_outer);
	std::vector<bool> is_point_within(n_points(), false);
	for (auto p: p_inner) 
		is_point_within[static_cast<GridPoint*>(p)->ind]=true;

	//2) copy cells info with only outer points
	//their point pointers temporary refer to this->points
	GridGeom res;
	for (int i=0; i<n_cells(); ++i){
		Cell c(res.n_cells());
		for (int j=0; j<cells[i]->dim(); ++j){
			int ind=cells[i]->points[j]->ind;
			if (!is_point_within[ind]){
				c.points.push_back(points[ind].get());
			}
		}
		if (c.dim()>2) aa::add_shared(res.cells, c);
	}

	//3) make a deep copy of all points which present in res.cells.points
	std::map<int, GridPoint*> inserted;
	for (int i=0; i<res.n_cells(); ++i){
		for (int j=0; j<res.cells[i]->dim(); ++j){
			GridPoint* this_pnt = res.cells[i]->points[j];
			auto fnd = inserted.find(this_pnt->ind);
			if (fnd!=inserted.end()){
				res.cells[i]->points[j]=fnd->second;
			} else {
				auto newp = aa::add_shared(res.points, *this_pnt);
				res.cells[i]->points[j]=newp;
				inserted.emplace(this_pnt->ind, newp);
			}
		}
	}

	//4) index nodes and return
	res.set_indicies();
	return res;
}

void GridGeom::force_cells_ordering(){
	for (auto c: cells) c->check_ordering();
}

shp_vector<GridGeom> GridGeom::subdivide() const{
	//cell->cell connectivity
	std::vector<std::vector<int>> cell_cell(n_cells());
	auto edges = get_edges();
	for (auto& e: edges){
		if (e.cell_left>=0 && e.cell_right>=0){
			cell_cell[e.cell_left].push_back(e.cell_right);
			cell_cell[e.cell_right].push_back(e.cell_left);
		}
	}
	//vector of cells use status (0 - unused, 1 - used)
	std::vector<int> cind(n_cells(), 0);

	//recursive algorithm of connected cells traversal
	std::function<void(int, vector<int>&)> add_cell = [&](int i, vector<int>& v ){
		cind[i]=1;
		v.push_back(i);
		for (auto ci: cell_cell[i]){
			if (cind[ci]==0) add_cell(ci, v);
		}
	};

	//build sets of single connected cells
	vector<vector<int>> sc_cells;
	while (1){
		auto fnd0 = std::find(cind.begin(), cind.end(), 0);
		//stop if all cells are already used
		if (fnd0 == cind.end()) break;
		sc_cells.push_back(vector<int>());
		add_cell(fnd0-cind.begin(), sc_cells.back());
	}

	//assemble new grids
	shp_vector<GridGeom> ret;
	for (auto& sc: sc_cells){
		ret.push_back(std::shared_ptr<GridGeom>(new GridGeom()));
		auto r = ret.back().get();
		r->add_data(*this, sc);
	}
	return ret;
}

GridGeom* GridGeom::combine(GridGeom* gmain, GridGeom* gsec){
	//1) build grids contours
	auto maincont = gmain->get_contours_collection();
	auto seccont = gsec->get_contours_collection();
	//2) input date to wireframe format
	PtsGraph wmain(*gmain);
	PtsGraph wsec(*gsec);
	//3) cut outer grid with inner grid contour
	wmain = PtsGraph::cut(wmain, seccont, -1);
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
	return ret;
}

GridGeom* GridGeom::cross_grids(GridGeom* gmain, GridGeom* gsec, 
		double buffer_size, double density, bool preserve_bp, crossgrid_callback cb){
	//initial scaling before doing anything
	if (cb("Building grid cross", "Scaling", 0, -1) == CALLBACK_CANCEL) return 0;
	auto sc = gmain->do_scale();
	gsec->do_scale(sc);
	buffer_size/=sc.L;
	
	//1 ---- combine grids without using buffers
	if (cb("Building grid cross", "Combining", 0.05, -1) == CALLBACK_CANCEL) return 0;
	std::unique_ptr<GridGeom> comb(GridGeom::combine(gmain, gsec));

	//2 ---- fill buffer zone
	//zero buffer requires no further actions
	if (buffer_size>geps){
	
		//loop over each get secondary contour
		auto csec  = gsec->get_contours();

		//loop over each single connected grid
		auto sg = comb->subdivide();

		//callback percentages 
		double cb_step = 1.0/(sg.size()*csec.size());
		double cb_cur = 0;
		for (auto grid: sg){
			for (auto c: csec){
				cb_cur+=cb_step;
				//1. filter out a grid from buffer zone for the contour
				if (cb("Building grid cross", "Filling buffer", 0.5, cb_cur-cb_step) 
						== CALLBACK_CANCEL) return 0;
				BufferGrid bg(*grid, c, buffer_size);

				//continue if buffergrid is empty
				if (bg.num_orig_cells()==0) continue; 

				//2. perform triangulation of buffer grid area
				if (cb("Building grid cross", "Filling buffer", 0.5, cb_cur-0.8*cb_step)
						== CALLBACK_CANCEL) return 0;
				auto bgcont = bg.boundary_info(preserve_bp);
				TriGrid g3(std::get<0>(bgcont), std::get<1>(bgcont), density);

				//3. change the internal of bg by g3ref grid
				if (cb("Building grid cross", "Filling buffer", 0.5, cb_cur-0.6*cb_step)
						== CALLBACK_CANCEL) return 0;
				bg.change_internal(g3);
				
				//4. update original grid using new filling of buffer grid
				if (cb("Building grid cross", "Filling buffer", 0.5, cb_cur-0.3*cb_step)
						== CALLBACK_CANCEL) return 0;
				bg.update_original();
			}
		}
		if (cb("Building grid cross", "Filling buffer", 0.5, 1.0)
						== CALLBACK_CANCEL) return 0;

		if (cb("Building grid cross", "Merging", 0.85, -1)
						== CALLBACK_CANCEL) return 0;
		//connect single connected grids
		comb->clear();
		for (size_t i=0; i<sg.size(); ++i){
			comb->add_data(*sg[i]);
		}
		if (sg.size()>0) comb->merge_congruent_points();
	}

	//scale back after all procedures have been done and return
	cb("Building grid cross", "Unscaling", 0.95, -1);
	gmain->undo_scale(sc);
	gsec->undo_scale(sc);
	comb->undo_scale(sc);

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
		for (int i=0; i<c->dim(); ++i){
			auto p1 = c->get_point(i-1);
			auto p2 = c->get_point(i);
			auto p3 = c->get_point(i+1);
			if (fabs(triarea(*p1, *p2, *p3))>geps){
				ret.push_back((*p1 + *p2 + *p3)/3.0);
				break;
			}
		}
	}
	return ret;
}




