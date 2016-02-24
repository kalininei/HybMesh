#include <map>
#include "addalgo.hpp"
#include "buffergrid.h"
namespace{

//this is a temporary function which will be deleted after
//crossgrid::Contour class will be completely substituted by HMCont2D classes
HMCont2D::Container<HMCont2D::Contour> contour_to_hm(const Contour& c){
	std::vector<Point> pnt;
	for (int i=0; i<c.n_points(); ++i) pnt.push_back(*c.get_point(i));
	return HMCont2D::Constructor::ContourFromPoints(pnt, true);
}
}

BufferGrid::BufferGrid(GridGeom& main, const PContour& cont, double buffer_size):GridGeom(){
	orig = &main;
	buffer = buffer_size;
	source_cont = Contour(cont);

	//check if cont and main intersection meets requirements for building the buffer
	auto treeorig = GGeom::Info::Contour(*orig);
	auto hmcont = contour_to_hm(source_cont);
	HMCont2D::Container<HMCont2D::ContourTree> cross;
	if (source_cont.area() > 0) cross = HMCont2D::Clip::Difference(treeorig, hmcont);
	else cross = HMCont2D::Clip::Intersection(treeorig, hmcont);
	HMCont2D::Clip::Heal(cross);
	if (HMCont2D::Area(cross)<geps) return;
	
	//1) find measures to all points
	vector<const Point*> pp;
	for (int i=0; i<main.n_points(); ++i) pp.push_back(main.get_point(i));
	auto meas = cont.meas_points(pp);

	//2) some other checks. I'm not sure why it is here. It was there before i've added intersection check.
	//May be these are artefacts from very old algos or smth. test07 falls without it.
	//Consider removing with refactoring.
	{
		//if any point the main grid lies strongly within the contour -> proceed with the procedure
		auto fnd1 = std::find_if(meas.begin(), meas.end(), [](double x){ return x>geps2; });
		if (fnd1 != meas.end()) goto PROCEED_PROC;
	}
	{
		//if less then two points of main lie on the cont -> stop procedure
		int bndcount = std::count_if(meas.begin(), meas.end(), [](double x){ return fabs(x)<geps2; });
		if (bndcount<2) return;
	}
	{
		//if any boundary main edge lies on cont -> proceed
		for (auto e: main.get_edges()) if (e.is_boundary()){
			if (fabs(meas[e.p1])<geps2 && fabs(meas[e.p2])<geps2) goto PROCEED_PROC;
		}
	}
	PROCEED_PROC:;

	//3) find all cells which includes filtered points
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

	//4) filter out all points which presents in inc_cells
	std::set<const Point*> points_set;
	for (auto c: orig_cells){
		for (int j=0; j<c->dim(); ++j) points_set.insert(c->get_point(j));
	}

	//5) build a grid from extracted cells
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

//boundary_info helper functions
namespace{

void simplify_bnd_edges(HMCont2D::ContourTree& cont, const HMCont2D::ECollection& bedges){
	std::vector<const Point*> not_needed_points;
	auto etree = HMCont2D::ExtendedTree::Assemble(bedges);
	for (int i=0; i<etree.cont_count(); ++i){
		auto ap = etree.get_contour(i)->all_points();
		auto cp = etree.get_contour(i)->corner_points();
		std::set<const Point*> cps(cp.begin(), cp.end());
		for (auto a: ap){
			if (cps.find(a) == cps.end()) not_needed_points.push_back(a);
		}
	}
	cont.RemovePoints(not_needed_points);
}

//extracts edges of 'from' which lie on 'where'
HMCont2D::ECollection extract_edges(const HMCont2D::ECollection& from,
		const HMCont2D::ECollection& where){
	std::map<const Point*, double> pdist;
	for(auto p: from.all_points()){
		pdist[p] = std::get<1>(HMCont2D::ECollection::FindClosestEdge(where, *p));
	}
	HMCont2D::ECollection ret;
	for(auto e: from){
		double dist1 = pdist[e->pstart];
		double dist2 = pdist[e->pend];
		if (ISZERO(dist1) && ISZERO(dist2)) ret.add_value(e);
	}
	return ret;
}

}//namespace

//1) get contour from buffer grid
//2) mark all edges as built on source contour/original grid contour/internal grid edges
//3) from resulting contour delete all non-node points from source contour
//4)                               non-corner points from original grid if preserve_bp==false
//5) for all points calculate weights as maximum distance to the left/right non-bnd points along contour
//6) for bnd points calculate weight using closest source and closest internal edge length
std::tuple<
	HMCont2D::ContourTree,
	std::map<Point*, double>
> BufferGrid::boundary_info(bool preserve_bp) const{
	//prepare return
	std::tuple<
		HMCont2D::ContourTree,
		std::map<Point*, double>
	> ret;
	auto& cont = std::get<0>(ret);
	auto& bw = std::get<1>(ret);

	//outer contour initial
	cont = GGeom::Info::Contour(*this);

	// ================ building edge types
	//edges types:
	//  1 - source contour edges, 
	//  2 - original grid boundary edges
	//  3 - original grid internal edges 
	auto scont = contour_to_hm(source_cont);
	auto grid_cont = GGeom::Info::Contour(*orig);
	HMCont2D::ECollection inner_edges = extract_edges(cont, scont);
	HMCont2D::ECollection bnd_edges = extract_edges(cont, grid_cont);
	HMCont2D::ECollection outer_edges;
	std::map<const HMCont2D::Edge*, int> edtypes;
	for (auto& e: cont) edtypes[e.get()] = 3;
	for (auto& e: inner_edges) edtypes[e.get()] = 1;
	for (auto& e: bnd_edges) edtypes[e.get()] = 2;
	
	//inner, bnd, outer edges set. Edges are shared with cont
	for (auto& ed: cont.data){
		if (edtypes[ed.get()] == 3) outer_edges.add_value(ed);
	}
	assert(inner_edges.size()!=0 || outer_edges.size()!=0);

	// ================ cont purge procedures
	//remove source contour points which doesn't equal source contour nodes
	//to get rid of hanging nodes
	std::vector<const Point*> bad_points;
	for (auto c: cont.nodes){
		auto cinfo = c->ordered_info(); cinfo.resize(cinfo.size()-1);
		for (auto& ci: cinfo){
			if (edtypes[ci.eprev.get()] == 1 && edtypes[ci.enext.get()] == 1){
				auto dst = HMCont2D::PCollection::FindClosestNode(scont.pdata, *ci.p);
				if (!ISZERO(sqrt(std::get<2>(dst)))) bad_points.push_back(ci.p);
			}
		}
	}
	cont.RemovePoints(bad_points);
	
	//simplify bnd edges if needed:
	//  for connected sequence of bnd_edges changes last point of first edge
	//  and deletes all others from contour tree
	if (!preserve_bp) simplify_bnd_edges(cont, bnd_edges);

	//Edges deleted from cont still present in bnd_edges collection but have NULL data.
	//clear *_edges ecollections from NULL data edges for further weight calculations
	auto rmfun = [](shared_ptr<HMCont2D::Edge> e){return e->pstart == 0;};
	inner_edges.data.resize(std::remove_if(inner_edges.begin(), inner_edges.end(), rmfun)-inner_edges.begin());
	outer_edges.data.resize(std::remove_if(outer_edges.begin(), outer_edges.end(), rmfun)-outer_edges.begin());
	bnd_edges.data.resize(std::remove_if(bnd_edges.begin(), bnd_edges.end(), rmfun)-bnd_edges.begin());

	// ========================= calculate weight procedures
	//set weights for all edges points (-1 for boundaries to process it further)
	auto calc_weight = [](double len_prev, double len_next, int tp_prev, int tp_next)->double{
			if (tp_prev == 1 && tp_next == 1) return std::max(len_prev, len_next);
			if (tp_prev == 1 && tp_next == 2) return len_prev;
			if (tp_prev == 1 && tp_next == 3) return std::max(len_prev, len_next);
			if (tp_prev == 2 && tp_next == 1) return len_next;
			if (tp_prev == 2 && tp_next == 2) return -1;
			if (tp_prev == 2 && tp_next == 3) return len_next;
			if (tp_prev == 3 && tp_next == 1) return std::max(len_prev, len_next);
			if (tp_prev == 3 && tp_next == 2) return len_prev;
			if (tp_prev == 3 && tp_next == 3) return std::max(len_prev, len_next);
			return -1;
		};
	for (auto& c: cont.nodes){
		auto oinfo = c->ordered_info();
		auto olens = HMCont2D::ECollection::ELengths(*c);
		for (int i=0; i<c->size(); ++i){
			double iprev = (i==0)?c->size()-1:i-1;
			const HMCont2D::Edge *eprev=oinfo[i].eprev.get();
			const HMCont2D::Edge *enext=oinfo[i].enext.get();
			bw[oinfo[i].p] = calc_weight(olens[iprev], olens[i], edtypes[eprev], edtypes[enext]);
		}
	}

	//calculate weights for points which have no distance set yet:
	//   boundary points of original grid
	for (auto& d: bw) if (d.second<0){
		auto d1 = HMCont2D::ECollection::FindClosestEdge(inner_edges, *d.first);
		auto d2 = HMCont2D::ECollection::FindClosestEdge(outer_edges, *d.first);
		if (std::get<0>(d1) == 0) d.second = std::get<0>(d2)->length();
		else if (std::get<0>(d2) == 0) d.second = std::get<0>(d1)->length();
		else{
			double dist1 = std::get<1>(d1), dist2 = std::get<1>(d2);
			double len1 = std::get<0>(d1)->length(), len2 = std::get<0>(d2)->length();
			double w1 = dist1/(dist1 + dist2);
			d.second = (1 - w1)*len1 + w1*len2;
		}
	}

	//return
	return ret;
}

void BufferGrid::update_original() const{
	//1) backup original contour
	auto cont_orig = PointsContoursCollection(orig->get_contours_collection());

	//2) delete all original cells
	std::set<int> _ds;
	aa::Cfill_container(orig_cells, _ds, [](const Cell* p){ return p->get_ind(); });
	aa::remove_entries(orig->cells, _ds);

	//3) backup original boundary points
	std::set<const GridPoint*> bpts_orig = orig->get_bnd_points();

	//4) add all buffer cells and nodes
	//to the end of original data arrays
	orig->add_data(*this);

	//5) get all new boundary points
	std::set<const GridPoint*> bpts_tmp = orig->get_bnd_points();
	std::set<const GridPoint*> bpts_new;
	for (auto x: bpts_tmp)
		if (bpts_orig.find(x)==bpts_orig.end())
			bpts_new.insert(x);

	//6) find all congruent contour points
	//original point index -> buffer point index
	std::map<int, int> cong_p;
	std::vector<const Point*> nocong;
	for (auto j=bpts_new.begin(); j!=bpts_new.end(); ++j){
		bool find_cong = false;
		for (auto i=bpts_orig.begin(); i!=bpts_orig.end(); ++i){
			if (**i == **j){
				cong_p.emplace((*i)->get_ind(), (*j)->get_ind());
				find_cong = true;
				break;
			}
		}
		if (!find_cong) nocong.push_back(*j);
	}
	
	//7) check: all points which have no congruent point should lie on the
	//original contour (not within)
	auto flt = cont_orig.filter_points_i(nocong);
	if (std::get<1>(flt).size() != nocong.size())
		throw std::runtime_error("Buffer grid does not match original");
	
	//8) remove points of the original grid which lie on the source contour.
	//   this is essential since some original points may lie on the source contour
	//   but not in the buffer grid.
	auto buf_cnt = get_contours();
	std::set<int> bad_pts;
	for (auto bp = bpts_orig.begin(); bp!=bpts_orig.end(); ++bp){
		//we don't care if contour is inner or outer. Let it always be inner
		bool hint = true;
		//if point is congruent to point in buffer grid then
		//it is on the edge of source contour.
		if (cong_p.find((*bp)->get_ind())==cong_p.end()){
			for (auto& c: buf_cnt){
				int r = c.is_inside(**bp, &hint);
				if (r==0) {
					bad_pts.insert((*bp)->get_ind());
					break;
				}
			}
		}
	}

	//9) merge congruent points and remove bad points
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
	if (bad_pts.size()>0){
		for (int ic=0; ic<orig->n_cells(); ++ic){
			auto cell = orig->cells[ic].get();
			for (int j=0; j<cell->dim(); ++j){
				int ind = cell->get_point(j)->get_ind();
				if (bad_pts.find(ind) != bad_pts.end()){
					delete_point_of_cell(cell, j--);
				}
			}
		}
	}

	//10) delete unused
	orig->delete_unused_points();
	orig->set_indicies();
}

