#include <map>
#include "addalgo.hpp"
#include "buffergrid.h"
#include "fileproc.h"

BufferGrid::BufferGrid(GridGeom& main, const PContour& cont, double buffer_size):GridGeom(){
	orig = &main;
	buffer = buffer_size;
	source_cont = Contour(cont);
	//1) find measures to all points
	vector<const Point*> pp;
	for (int i=0; i<main.n_points(); ++i) pp.push_back(main.get_point(i));
	auto meas = cont.meas_points(pp);

	//2) find all cells which includes filtered points
	for (int i=0; i<main.n_cells(); ++i){
		auto c = main.get_cell(i);
		bool is_good = false;
		//if any point lies inside contour ignore cell
		//if any point lies inside outer buffer zone but not on contour add cell
		for (int j=0; j<c->dim(); ++j){
			double& m = meas[c->get_point(j)->get_ind()];
			//point on contour are ambiguous
			if (fabs(m)<geps2){
				is_good = true; continue;
			//if point in buffer zone add cell
			} else if (-m>0 && -m<sqr(buffer_size-geps) ){
				is_good=true; break;
			//if point within contour ignore cell
			} else if (m>0){
			       is_good=false; break; 
			}
		}
		if (is_good) orig_cells.push_back(c);
	}

	//3) filter out all points which presents in inc_cells
	std::set<const Point*> points_set;
	for (auto c: orig_cells){
		for (int j=0; j<c->dim(); ++j) points_set.insert(c->get_point(j));
	}
	
	//4) build a grid from extracted cells
	//points
	std::map<int, GridPoint*> mp;
	for(auto p: points_set){
		aa::add_shared(points, GridPoint(*p));
		mp[static_cast<const GridPoint*>(p)->get_ind()]=points.back().get();
	}
	//cells
	for (auto c: orig_cells){
		auto newc = aa::add_shared(cells, Cell());
		for (int i=0; i<c->dim(); ++i){
			add_point_to_cell(newc, mp[c->get_point(i)->get_ind()]);
		}
	}
	set_indicies();
}

std::tuple<
	std::set<const Point*>,   //inner_bp - boundary points from inner contour
	std::set<const Point*>,   //outer_bp - boundary points from outer contour
	std::set<const Point*>    //true_bp  - boundary points from original boundary
> BufferGrid::build_bedges() const{
	//return value initializatin
	std::tuple<
		std::set<const Point*>,
		std::set<const Point*>,
		std::set<const Point*>
	> ret;
	auto& inner_bp = std::get<0>(ret);
	auto& outer_bp = std::get<1>(ret);
	auto& true_bp = std::get<2>(ret);

	//get contour points
	auto cnts = get_contours();
	vector<const Point*> pp;
	for (auto c: cnts) for (int i=0; i<c.n_points(); ++i){
		pp.push_back(c.get_point(i));
	}
	//get measures of contour points
	auto meas = source_cont.meas_points(pp);
	//sort points depending on their meassures
	for (size_t i=0; i<meas.size(); ++i){
		//if point lies without buffer zone this is outer contour
		if (-meas[i]>sqr(buffer-geps)) outer_bp.insert(pp[i]);
		//if point lies on the contour this is inner contour point
		else if (fabs(meas[i])<geps2) inner_bp.insert(pp[i]);
		//if this is boundary point from the original grid
		else true_bp.insert(pp[i]);
	}
	return ret;
}

void BufferGrid::new_edge_points(Edge& e, const vector<double>& wht){
	auto p1 = points[e.p1].get(), p2 = points[e.p2].get();
	vector<GridPoint*> newp;
	for (auto w: wht){
		newp.push_back(aa::add_shared(points, GridPoint(Point::Weigh(*p1, *p2, w))));
	}
	//add to cell_left
	if (e.cell_left>=0){
		auto cleft = cells[e.cell_left].get();
		Cell* new_cleft = new Cell();
		for (int i=0; i<cleft->dim(); ++i){
			auto cell_point = cleft->get_point(i);
			add_point_to_cell(new_cleft, const_cast<GridPoint*>(cell_point));
			if (cell_point==p1){
				for (auto p: newp) add_point_to_cell(new_cleft, p);
			}
		}
		cells[e.cell_left] = std::shared_ptr<Cell>(new_cleft);
	}
	//add to cell_right
	if (e.cell_right>=0){
		std::reverse(newp.begin(), newp.end());
		auto cright = cells[e.cell_right].get();
		Cell* new_cright = new Cell();
		for (int i=0; i<cright->dim(); ++i){
			auto cell_point = cright->get_point(i);
			add_point_to_cell(new_cright, const_cast<GridPoint*>(cell_point));
			if (cell_point==p2){
				for (auto p: newp) add_point_to_cell(new_cright, p);
			}
		}
		cells[e.cell_right] = std::shared_ptr<Cell>(new_cright);
	}
}

std::tuple<
	std::vector<PContour>,
	std::vector<double>
> BufferGrid::boundary_info() const{
	//return value
	std::tuple<vector<PContour>, vector<double>> ret;
	auto& real_contour = std::get<0>(ret);
	auto& lc = std::get<1>(ret);

	//nodes origin info
	auto origin = build_bedges();
	auto& inner_bp = std::get<0>(origin);
	auto& outer_bp = std::get<1>(origin);
	auto& true_bp = std::get<2>(origin);
	// --- build real contour
	real_contour = get_contours();
	//delete points which lie on the sections of the source contour
	for (auto& c: real_contour){
		std::set<int> bad_points;
		for (int i=0; i<c.n_points(); ++i){
			if (!c.is_corner_point(i) && inner_bp.find(c.get_point(i))!=inner_bp.end()){
				auto p = c.get_point(i);
				//find if p is the edge of source contour
				for (int j=0; j<source_cont.n_points(); ++j){
					if (*p == *source_cont.get_point(j)){
						p=0; break;
					}
				}
				if (p!=0) bad_points.insert(i);
			}
		}
		c.delete_by_index(bad_points);
	}

	//build characteristic steps
	int _npt=0; for (auto& c: real_contour) _npt+=c.n_points();
	lc.reserve(_npt);  //dist_dict points to lc data hence reserve is necessary
	std::map<const Point*, double*> dist_dict;
	for (auto& c: real_contour){
		//types of vertices: inner, outer, true boundary
		vector<int> types;
		for (int i=0; i<c.n_points(); ++i){
			if (inner_bp.find(c.get_point(i))!=inner_bp.end()) types.push_back(0);
			else if (true_bp.find(c.get_point(i))!=true_bp.end()) types.push_back(2);
			else types.push_back(1);
		}
		//calculating distances
		vector<double> dist = c.section_lenghts();
		//characteristic steps
		for (int i=0; i<c.n_points(); ++i){
			int iprev = (i==0)?c.n_points()-1:i-1;
			int inext = (i==c.n_points()-1)?0:i+1;
			double d = (dist[iprev]+dist[i])/2.0 ;
			if (types[i]==types[inext] && types[i]!=types[iprev]) d = dist[i];
			else if (types[i]!=types[inext] && types[i]==types[iprev]) d = dist[iprev];
				
			lc.push_back(d);
			dist_dict.emplace(c.get_point(i), &lc.back());
		}
	}
	
	//for true boundary points find closest inner and outer points
	auto closest_point = [](const Point* pnt, const std::set<const Point*>& col)
			->std::pair<const Point*, double>{
		double dist = gbig;
		const Point* closest = 0;
		for (auto& p2: col){
			double d = Point::meas(*pnt, *p2);
			if (d<dist) { dist = d; closest = p2; }
		}
		return std::make_pair(closest, sqrt(dist));
	};
	for (auto& p: true_bp){
		auto closest_inner = closest_point(p, inner_bp);
		auto closest_outer = closest_point(p, outer_bp);
		double f1 = *dist_dict[closest_inner.first], f2 = *dist_dict[closest_outer.first];
		double d1 = closest_inner.second, d2 = closest_outer.second;
		*dist_dict[p] = (f1*d1+f2*d2)/(d1+d2);
	}

	return ret;
}

void BufferGrid::update_original() const{
	//1) delete all original cells
	std::set<int> _ds;
	aa::Cfill_container(orig_cells, _ds, [](const Cell* p){ return p->get_ind(); });
	aa::remove_entries(orig->cells, _ds);
	//2) add all buffer cells and nodes
	orig->add_data(*this);
	//3) find all congruent contour points
	orig->set_indicies();
	std::map<int, int> cong_p;
	auto ocont = orig->get_contours();
	std::vector<const GridPoint*> bpts;
	for (auto c: ocont){
		for (int i=0; i<c.n_points(); ++i) bpts.push_back(static_cast<const GridPoint*>(c.get_point(i)));
	}
	for (int i=0; i<bpts.size(); ++i){
		for (int j=i+1; j<bpts.size(); ++j){
			if (*bpts[i] == *bpts[j]){
				cong_p.emplace(bpts[i]->get_ind(), bpts[j]->get_ind());
				break;
			}
		}
	}
	//4) merge congruent points
	for (int ic=0; ic<orig->n_cells(); ++ic){
		auto cell = orig->cells[ic].get();
		for (int j=0; j<cell->dim(); ++j){
			int ind = cell->get_point(j)->get_ind();
			auto fnd = cong_p.find(ind);
			if (fnd!=cong_p.end()){
				change_point_of_cell(cell, j, orig->points[fnd->second].get());
			}
		}
	}
	//5) delete unused 
	orig->delete_unused_points();
	orig->set_indicies();
}

