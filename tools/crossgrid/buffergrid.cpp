#include <map>
#include "addalgo.hpp"
#include "buffergrid.h"

BufferGrid::BufferGrid(GridGeom& main, const PContour& cont, double buffer_size):GridGeom(){
	orig = &main;
	buffer_size*=buffer_size; //quadratic comparison is used
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
			} else if (m<0 && m>-buffer_size+geps2){
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
	std::map<int, GPoint*> mp;
	for(auto p: points_set){
		aa::add_shared(points, GPoint(*p));
		mp[static_cast<const GPoint*>(p)->get_ind()]=points.back().get();
	}
	//cells
	for (auto c: orig_cells){
		aa::add_shared(cells, Cell());
		for (int i=0; i<c->dim(); ++i){
			add_point_to_cell(cells.back().get(), mp[c->get_point(i)->get_ind()]);
		}
	}
	set_indicies();
}

void BufferGrid::split_bedges(){
	//TODO
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
	std::vector<const GPoint*> bpts;
	for (auto c: ocont){
		for (int i=0; i<c.n_points(); ++i) bpts.push_back(static_cast<const GPoint*>(c.get_point(i)));
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



