#include <assert.h>
#include <map>
#include "grid.h"
#include "addalgo.hpp"
#include "fileproc.h"
#include "trigrid.h"
#include "buffergrid.h"
#include "wireframegrid.h"

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

GridGeom* GridGeom::combine(GridGeom* gmain, GridGeom* gsec){
	//1) input date to wireframe format
	PtsGraph wmain(*gmain);
	PtsGraph wsec(*gsec);
	//2) cut outer grid with inner grid contour
	wmain = PtsGraph::cut(wmain, gsec->get_contours(),-1);
	//3) overlay grids
	wmain = PtsGraph::overlay(wmain, wsec);
	//4) return
	return new GridGeom(wmain.togrid());
}

GridGeom* GridGeom::cross_grids(GridGeom* gmain, GridGeom* gsec, double buffer_size, double density){
	//initial scaling before doing anything
	auto sc = gmain->do_scale();
	gsec->do_scale(sc);
	buffer_size/=sc.L;
	
	//1 ---- combine grids without using buffers
	GridGeom* comb = GridGeom::combine(gmain, gsec);

	//2 ---- fill buffer zone
	//zero buffer requires no further actions
	if (buffer_size>geps){
		//loop over each get secondary contour
		auto csec  = gsec->get_contours();

		for (auto c: csec){
			//1. filter out a grid from buffer zone for the contour
			BufferGrid bg(*comb, c, buffer_size);

			//2. perform triangulation of buffer grid area
			auto bgcont = bg.boundary_info();
			TriGrid g3(std::get<0>(bgcont), std::get<1>(bgcont), density);

			//3. change the internal of bg by g3ref grid
			bg.change_internal(g3);
			
			//4. update original grid using new filling of buffer grid
			bg.update_original();
		}
	}

	//scale back after all procedures have been done and return
	gmain->undo_scale(sc);
	gsec->undo_scale(sc);
	comb->undo_scale(sc);
	return comb;
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
