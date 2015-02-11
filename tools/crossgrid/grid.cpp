#include <assert.h>
#include <map>
#include "grid.h"
#include "addalgo.hpp"
#include "fileproc.h"
#include "trigrid.h"
#include "buffergrid.h"

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


void GridGeom::add_data(const GridGeom& g){
	int start_p_index = n_points();
	for (auto p: g.points) aa::add_shared(points, GPoint(*p));
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
		aa::add_shared(points, GPoint(x,y,i));
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
		is_point_within[static_cast<GPoint*>(p)->ind]=true;

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
	std::map<int, GPoint*> inserted;
	for (int i=0; i<res.n_cells(); ++i){
		for (int j=0; j<res.cells[i]->dim(); ++j){
			GPoint* this_pnt = res.cells[i]->points[j];
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

GridGeom* GridGeom::cross_grids(GridGeom* gmain, GridGeom* gsec, double buffer_size){
	//normalize
	auto nrm = gmain->do_scale();
	gsec->do_scale(nrm);

	//get all contours from secondary grid
	auto csec  = gsec->get_contours();

	//delete each contour area from the main grid
	GridGeom constr(*gmain);
	for (auto& c2: csec){
		auto cont = c2.widen_contour(buffer_size);
		constr=constr.remove_area(cont);
	}

	//after grid was cleared from area wich contains secondary grid
	//we can add it
	constr.add_data(*gsec);

	//renormalize
	gmain->undo_scale(nrm);
	gsec->undo_scale(nrm);
	constr.undo_scale(nrm);

	//return grid sum
	return new GridGeom(constr);
}

GridGeom* GridGeom::combine(GridGeom* gmain, GridGeom* gsec){
	//TODO
	//------ DUMMY output: [0,0]x[1,1] + [0.3, 0.3]x[0.6, 0.6]
	GridGeom* ret = new GridGeom(*gsec);
	for (int i=0; i<gmain->n_points(); ++i){
		aa::add_shared(ret->points, GPoint(*gmain->get_point(i)));
	}
	for (int i=0; i<gmain->n_cells(); ++i){
		auto c = gmain->get_cell(i);
		auto cnew = Cell();
		for (int j=0; j<c->dim(); ++j){
			int ind = c->get_point(j)->ind;
			cnew.points.push_back(ret->points[ind+gsec->n_points()].get());
		}
		aa::add_shared(ret->cells, Cell(cnew));
	}
	int bc = gsec->n_cells();
	//--- corner cells
	//left bottom cell
	ret->cells[bc+22]->points[2] = ret->points[0].get();
	//right bottom cell
	ret->cells[bc+26]->points[3] = ret->points[30].get();
	//top left
	ret->cells[bc+62]->points[1] = ret->points[31*30].get();
	//top right
	ret->cells[bc+66]->points[0] = ret->points[31*31-1].get();
	////--- left cells
	for (int i=0; i<3; ++i){
		int cind = bc+32+10*i;
		for (int j=0; j<11; ++j){
			ret->cells[cind]->points.push_back(ret->points[31*(10*i +j)].get());
		}
		ret->cells[cind]->points.push_back(ret->cells[cind]->points[3]);
		aa::remove_entries(ret->cells[cind]->points, {1,2,3});
	}
	//bottom cells
	for (int i=0; i<3; ++i){
		int cind = bc+23+i;
		for (int j=0; j<11; ++j){
			ret->cells[cind]->points.push_back(ret->points[10*(i+1)-j].get());
		}
		aa::remove_entries(ret->cells[cind]->points, {2,3});
	}
	//right cells
	for (int i=0; i<3; ++i){
		int cind = bc+36+10*i;
		for (int j=0; j<11; ++j){
			ret->cells[cind]->points.push_back(ret->points[ (10*(i+1)-j)*31+30 ].get());
		}
		aa::remove_entries(ret->cells[cind]->points, {0,3});
	}
	//top cells
	for (int i=0; i<3; ++i){
		int cind = bc+63+i;
		for (int j=0; j<11; ++j){
			ret->cells[cind]->points.push_back(ret->points[30*31+10*i+j].get());
		}
		ret->cells[cind]->points.push_back(ret->cells[cind]->points[2]);
		ret->cells[cind]->points.push_back(ret->cells[cind]->points[3]);
		aa::remove_entries(ret->cells[cind]->points, {0,1,2,3});
	}
	//superfluous cells
	std::set<int> sup_cells;
	for (int i = 3; i<6; ++i){
		for (int j=3; j<6; ++j){
			sup_cells.insert(bc+j*10+i);
		}
	}
	aa::remove_entries(ret->cells, sup_cells);
	ret->delete_unused_points();
	ret->set_indicies();
	std::cout<<"DUMMY COMBINE GRID: "<<ret->n_points()<<"  "<<ret->n_cells()<<std::endl;
	return ret;
}

GridGeom* GridGeom::cross_grids2(GridGeom* gmain, GridGeom* gsec, double buffer_size){
	//1 ---- combine grids without using buffers
	GridGeom* comb = GridGeom::combine(gmain, gsec);

	//2 ---- fill buffer
	//loop over each get secondary contour
	auto csec  = gsec->get_contours();

	for (auto c: csec){
		//1. filter out a grid from buffer zone for the contour
		BufferGrid bg(*comb, c, buffer_size);
		
		//2. add nodes to edges which are boundary in comb.first grid
		bg.split_bedges();

		//3. triangulate bg area
		std::vector<PContour> bg_cont = bg.get_contours();
		TriGrid  g3(bg_cont);

		//4 build a grid within the area
		TriGrid g3ref = g3.refine_grid(0.5);
	
		//5. change the internal of bg by g3ref grid
		bg.change_internal(g3ref);
		
		//6. update original grid using new filling of buffer grid
		bg.update_original();
	}

	return comb;
}

void GridGeom::change_internal(const GridGeom& gg){
	//----------- subprocedures
	auto boundary_points = [](const vector<PContour>& cv)->std::set<const GPoint*>{
		std::set<const GPoint*> ret;
		for (auto c: cv){
			for (int i=0; i<c.n_points(); 
				ret.insert(static_cast<const GPoint*>(c.get_point(i++))));
		}
		return ret;
	};
	auto internal_points = [](const GridGeom& grid, const std::set<const GPoint*>& bs)
			->std::set<const GPoint*>{
		std::set<const GPoint*> ret;
		for (int i=0; i<grid.n_points(); ++i){
			if (bs.find(grid.get_point(i))==bs.end())
				ret.insert(static_cast<const GPoint*>(grid.get_point(i)));

		};
		return ret;
	};
	auto boundary_match = [](const std::set<const GPoint*>& s1, std::set<const GPoint*>& s2)
			->std::map<const GPoint*, const GPoint*>{
		std::map<const GPoint*, const GPoint*> ret;
		for (auto p: s1){
			for (auto sit = s2.begin(); sit!=s2.end(); ++sit){
				if (*p==**sit){
					ret.emplace(p, *sit);
					s2.erase(sit);
					break;
				}
			}
		}
		assert(ret.size()==s1.size() && s2.size()==0);
		return ret;
	};
	auto fill_points = [](const GridGeom& g, const shp_vector<GPoint>& bs, 
			const std::set<const GPoint*>& is, std::map<const GPoint*, const GPoint*>& m)
			->shp_vector<GPoint>{
		shp_vector<GPoint> ret = bs;
		for (auto p: is) m[p] = aa::add_shared(ret, GPoint(*p));
		return ret;
	};
	auto fill_cells = [](const shp_vector<GPoint>& p, const std::map<const GPoint*, const GPoint*>& m,
			const GridGeom& grid)->shp_vector<Cell>{
		shp_vector<Cell> ret;
		for (int ic=0; ic<grid.n_cells(); ++ic){
			auto c = grid.get_cell(ic);
			auto newc = aa::add_shared(ret, Cell());
			for (int i=0; i<c->dim(); ++i){
				auto p = m.find(c->get_point(i))->second;
				add_point_to_cell(newc, const_cast<GPoint*>(p));
			}
		}
		return ret;
	};
	//------------ main procedure
	//1) build outer contours 
	auto cont1 = get_contours();
	auto cont2 = gg.get_contours();
	//2) create sets and vectors of boundary and internal points
	std::set<const GPoint*> bs1 = boundary_points(cont1);
	std::set<const GPoint*> bs2 = boundary_points(cont2);
	std::set<const GPoint*> is2 = internal_points(gg, bs2);
	shp_vector<GPoint> boundary = 
		aa::Cfill_vector(bs1, [this](const GPoint* p){ return this->points[p->get_ind()]; });
	//3) match between boundary sets
	std::map<const GPoint*, const GPoint*> s1s2 = boundary_match(bs2, bs1);
	//4) fill points array and add them to s1s2 dictionary
	shp_vector<GPoint> new_points = fill_points(*this, boundary, is2, s1s2);
	//5) fill cells array
	shp_vector<Cell> new_cells = fill_cells(new_points, s1s2, gg);
	//6) swap
	std::swap(points, new_points);
	std::swap(cells, new_cells);
	set_indicies();
}
