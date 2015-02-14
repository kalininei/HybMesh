#include <map>
#include "addalgo.hpp"
#include "buffergrid.h"
#include "fileproc.h"

BufferGrid::BufferGrid(GridGeom& main, const PContour& cont, double buffer_size):GridGeom(){
	orig = &main;
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
	std::map<int, GPoint*> mp;
	for(auto p: points_set){
		aa::add_shared(points, GPoint(*p));
		mp[static_cast<const GPoint*>(p)->get_ind()]=points.back().get();
	}
	//cells
	for (auto c: orig_cells){
		auto newc = aa::add_shared(cells, Cell());
		for (int i=0; i<c->dim(); ++i){
			add_point_to_cell(newc, mp[c->get_point(i)->get_ind()]);
		}
	}
	set_indicies();
	build_bedges(cont, buffer_size);
}

BufferGrid::BndEdge::BndEdge(const Edge& e): Edge(e.p1, e.p2){
	cell_left = e.cell_left; cell_right = e.cell_right;
	w1=-1; w2=-1;
}
void BufferGrid::build_bedges(const PContour& cont, double buffer_size){
	//inner and outer boundary points
	//boundary point -> boundary step dictionary
	std::map<const Point*, double> inner_bp;  //boundary points from inner contour
	std::map<const Point*, double> outer_bp;  //boundary points from outer contour
	std::map<const Point*, double> true_bp;   //boundary points from original boundary
	
	//get contour points
	auto cnts = get_contours();
	vector<const Point*> pp;
	for (auto c: cnts) for (int i=0; i<c.n_points(); ++i){
		pp.push_back(c.get_point(i));
	}
	//get measures of contour points
	auto meas = cont.meas_points(pp);
	//sort points depending on their meassures
	for (size_t i=0; i<meas.size(); ++i){
		//if point lies without buffer zone this is outer contour
		if (-meas[i]>sqr(buffer_size-geps)) outer_bp.emplace(pp[i], -1);
		//if point lies on the contour this is inner contour point
		else if (fabs(meas[i])<geps2) inner_bp.emplace(pp[i], -1);
		//if this is boundary point from the original grid
		else true_bp.emplace(pp[i], -1);
	}
	//fill steps for inner and outer points
	std::map<const Point*, double*> all_bp;
	for (auto& v: inner_bp) all_bp.emplace(v.first, &v.second);
	for (auto& v: outer_bp) all_bp.emplace(v.first, &v.second);
	for (auto& c: cnts){
		auto pprev = c.get_point(c.n_points()-1);
		for (int i=0; i<c.n_points(); ++i){
			auto p = c.get_point(i);
			auto fndp = all_bp.find(p);
			auto fndprev = all_bp.find(pprev);
			if (fndp!=all_bp.end() && fndprev!=all_bp.end()){
				double Len = Point::dist(*p, *pprev);
				if (*fndp->second>0) *fndp->second = (*fndp->second + Len)/2.0;
				else *fndp->second = Len;
				if (*fndprev->second>0) *fndprev->second = (*fndprev->second + Len)/2.0;
				else *fndprev->second = Len;
			}
			pprev = p;
		}
	}

	//fill steps for true_bp as a weighted value 
	//between steps of closest inner and closest outer point
	auto closest_point = [](const Point* pnt, const std::map<const Point*, double>& col)
			->std::pair<const Point*, double>{
		double dist = gbig;
		const Point* closest = 0;
		for (auto& p2: col){
			double d = Point::meas(*pnt, *p2.first);
			if (d<dist) { dist = d; closest = p2.first; }
		}
		return std::make_pair(closest, sqrt(dist));
	};
	for (auto& bp: true_bp){
		//find closest inner and outer points
		auto closest_inner = closest_point(bp.first, inner_bp);
		auto closest_outer = closest_point(bp.first, outer_bp);
		//find steps
		double f1 = inner_bp[closest_inner.first];
		double f2 = outer_bp[closest_outer.first];
		bp.second = (f1*closest_inner.second + f2*closest_outer.second)/(closest_inner.second + closest_outer.second);
	}
	//build edges
	for (auto& v: true_bp) all_bp.emplace(v.first, &v.second);
	auto edges = get_edges();
	for (auto& e: edges) if (e.is_boundary()){
		auto p1 = points[e.p1].get(), p2 = points[e.p2].get();
		//if any of p1 and p2 are in true_bp or
		//p1 and p2 are in different sets then ignore this edge
		bool p1_inner = (inner_bp.find(p1)!=inner_bp.end());
		bool p2_inner = (inner_bp.find(p2)!=inner_bp.end());
		if (p1_inner && p2_inner) continue;
		bool p1_outer = (outer_bp.find(p1)!=outer_bp.end());
		bool p2_outer = (outer_bp.find(p2)!=outer_bp.end());
		if (p1_outer && p2_outer) continue;
		//add edge
		bedges.push_back(BndEdge(e));
		bedges.back().w1 = *all_bp[p1];
		bedges.back().w2 = *all_bp[p2];
	}
}

void BufferGrid::new_edge_points(Edge& e, const vector<double>& wht){
	auto p1 = points[e.p1].get(), p2 = points[e.p2].get();
	vector<GPoint*> newp;
	for (auto w: wht){
		newp.push_back(aa::add_shared(points, GPoint(Point::Weigh(*p1, *p2, w))));
	}
	//add to cell_left
	if (e.cell_left>=0){
		auto cleft = cells[e.cell_left].get();
		Cell* new_cleft = new Cell();
		for (int i=0; i<cleft->dim(); ++i){
			auto cell_point = cleft->get_point(i);
			add_point_to_cell(new_cleft, const_cast<GPoint*>(cell_point));
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
			add_point_to_cell(new_cright, const_cast<GPoint*>(cell_point));
			if (cell_point==p2){
				for (auto p: newp) add_point_to_cell(new_cright, p);
			}
		}
		cells[e.cell_right] = std::shared_ptr<Cell>(new_cright);
	}
}

void BufferGrid::split_bedges(double density){
	for (auto e: bedges){
		auto p1 = *points[e.p1], p2 = *points[e.p2];
		double Len = Point::dist(p1, p2);
		auto refPnt = RefineSection(e.w1, e.w2, Len, density);
		aa::Cforeach(refPnt, [Len](double& x){x/=Len;});
		new_edge_points(e, refPnt);
	}
	set_indicies();
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



