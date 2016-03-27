#include <map>
#include "addalgo.hpp"
#include "buffergrid.h"
#include "procgrid.h"
#include "wireframegrid.h"
namespace{

//this is a temporary function which will be deleted after
//crossgrid::Contour class will be completely substituted by HMCont2D classes
HMCont2D::Container<HMCont2D::Contour> contour_to_hm(const Contour& c){
	std::vector<Point> pnt;
	for (int i=0; i<c.n_points(); ++i) pnt.push_back(*c.get_point(i));
	return HMCont2D::Constructor::ContourFromPoints(pnt, true);
}

HMCont2D::Container<HMCont2D::ECollection> get_contact(HMCont2D::ContourTree& gridcont, HMCont2D::Contour& source){
	PtsGraph g(gridcont), s(source);
	PtsGraph over = PtsGraph::overlay(g, s);
	//leave only those lines which centers lie on source
	//and inside gridcont
	vector<Point> cp = over.center_line_points();
	vector<int> srt = HMCont2D::Algos::SortOutPoints(gridcont, cp);
	std::set<int> bad_lines;
	for (int i=0; i<cp.size(); ++i) if (srt[i] != INSIDE) bad_lines.insert(i);
	over.exclude_lines(bad_lines);
	auto col = over.toecollection();
	return HMCont2D::Container<HMCont2D::ECollection>::DeepCopy(col);
}

//calculates measures to contact lines for points in pt.
//if points lies without bounding box of given contours + buf sets 2*buf*buf as its measure
vector<double> build_measures(HMCont2D::ECollection& src, const ShpVector<GridPoint>& pt){
	vector<double> ret;
	for (auto p: pt){
		auto r = HMCont2D::ECollection::FindClosestEdge(src, *p);
		double& av = std::get<1>(r);
		ret.push_back(av*av);
	}
	return ret;
}

//for points wich lie outside given area multiplies its measure by -1
void negative_to_outside(vector<double>& meas, const HMCont2D::Contour& area, const ShpVector<GridPoint>& pt){
	vector<Point> p2;
	for (auto p: pt) p2.push_back(*p);
	vector<int> srt = HMCont2D::Algos::SortOutPoints(area, p2);
	for (int i=0; i<meas.size(); ++i){
		if (srt[i] == OUTSIDE) meas[i]*=-1;
	}
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
	double domarea = HMCont2D::Area(cross);
	if (domarea < geps) return;

	//1) find measures to all points
	auto contact = get_contact(treeorig, hmcont);
	auto ptlist = GGeom::Info::SharePoints(main);
	auto meas = build_measures(contact, ptlist);
	negative_to_outside(meas, hmcont, ptlist);

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
				if (j!=c->dim()-1){ is_good = true; continue; }
				else {
					//all points lie on source contour. we have to check intersection
					auto ccont = GGeom::Info::CellContour(main, c->get_ind());
					double unarea = HMCont2D::Area(HMCont2D::Clip::Union(cross, ccont));
					if (fabs(unarea-domarea)>geps*geps) is_good = false;
					else is_good = true;
				}
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
	//std::set<const GridPoint*> points_set;
	//for (auto c: orig_cells){
		//for (int j=0; j<c->dim(); ++j) points_set.insert(c->get_point(j));
	//}
	std::set<int> points_set_ind;
	for (auto c: orig_cells){
		for (int j=0; j<c->dim(); ++j) points_set_ind.insert(c->get_point(j)->get_ind());
	}
	std::vector<const GridPoint*> points_set;
	for (auto ip: points_set_ind) points_set.push_back(main.get_point(ip));

	//5) build a grid from extracted cells
	//points
	std::map<int, GridPoint*> mp;
	for(auto p: points_set){
		aa::add_shared(points, GridPoint(*p));
		mp[p->get_ind()]=points.back().get();
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

//void simplify_bnd_edges(HMCont2D::ContourTree& cont,
//                std::map<const HMCont2D::Edge*, int>& edtypes, double angle0){
//        if (ISZERO(angle0)) angle0 = 1e-6;
//        //assemble needed points because there could be cases when a point
//        //is needed from one contour of 'cont' and not needed from the other
//        std::set<const Point*> needed_points;
//        for (auto c: cont.nodes){
//                vector<HMCont2D::Contour::PInfo> info = c->ordered_info();
//                for (auto i: info){
//                        bool fnd1 = edtypes[i.eprev.get()] == 2;
//                        bool fnd2 = edtypes[i.enext.get()] == 2;
//                        //point is not needed only if two adjacent edges are boundary and
//                        //form straight angle
//                        if (!fnd1 || !fnd2) needed_points.insert(i.p);
//                        else {
//                                //double area3 = triarea(*i.pprev, *i.p, *i.pnext);
//                                //if (fabs(area3)>geps*geps) needed_points.insert(i.p);
//                                double a = Angle(*i.pprev, *i.p, *i.pnext)/M_PI*180;
//                                if (fabs(a-180)>angle0) needed_points.insert(i.p);
//                        }
//                }
//        }
//        vector<const Point*> not_needed_points;
//        for (auto p: cont.all_points()){
//                if (needed_points.find(p) == needed_points.end()) not_needed_points.push_back(p);
//        }
//        cont.RemovePoints(not_needed_points);
//}

std::set<Point*> find_needed_points(HMCont2D::ContourTree& cont, HMCont2D::ECollection& bnded, double a0){
	a0 *= (M_PI/180);
	std::set<Point*> ret;
	for (auto c: cont.nodes){
		vector<HMCont2D::Contour::PInfo> info = c->ordered_info();
		for (auto i: info){
			if (!i.eprev || !i.enext ||
			    !bnded.contains(i.eprev.get()) || !bnded.contains(i.enext.get()) ||
			    fabs(Angle(*i.pprev, *i.p, *i.pnext)-M_PI)>a0){
				ret.insert(i.p);
			}
		}
	}
	return ret;
}

HMCont2D::ExtendedTree insignificant_segments(HMCont2D::ContourTree& cont, std::set<Point*>& needp){
	std::set<HMCont2D::Edge*> inserted_edges;
	HMCont2D::ECollection ret;
	for (auto c: cont.nodes){
		vector<HMCont2D::Contour::PInfo> info = c->ordered_info();
		for (auto i: info) if (needp.find(i.p) == needp.end()){
			auto insert1 = inserted_edges.emplace(i.eprev.get());
			auto insert2 = inserted_edges.emplace(i.enext.get());
			if (insert1.second) ret.add_value(i.eprev);
			if (insert2.second) ret.add_value(i.enext);
		}
	}
	auto rtree = HMCont2D::Assembler::ETree(ret);
	//match directions of all contours in rtree with directions in cont
	//!!! This could change directions of closed contours of
	//    ExtendedTree which are now directed according to its nested structure.
	//    However, we hope that nested structure of ExtendedTree
	//    will not be used later in the algorithm and this will not matter.
	for (auto c: rtree.all_contours()){
		assert(c->size() > 1);  //since we compose contours from point neighbors
		Point* p1 = c->first();
		Point* p2 = c->next_point(p1);
		HMCont2D::Contour* parcont = cont.get_contour(p1);
		Point* p3 = parcont->next_point(p1);
		if (p3 != p2) c->Reverse();
	}
	return rtree;
}
//HMCont2D::ExtendedTree insignificant_segments(HMCont2D::ContourTree& cont, HMCont2D::ECollection& bnded, double a0){
//        std::set<HMCont2D::Edge*> inserted_edges;
//        HMCont2D::ECollection ret;
//        for (auto c: cont.nodes){
//                vector<HMCont2D::Contour::PInfo> info = c->ordered_info();
//                for (auto i: info){
//                        if (!i.eprev || !i.enext) continue;
//                        if (!bnded.contains(i.eprev.get()) || !bnded.contains(i.enext.get())) continue;
//                        double a = Angle(*i.pprev, *i.p, *i.pnext)/M_PI*180;
//                        if (fabs(a-180)>a0) continue;
//                        auto insert1 = inserted_edges.emplace(i.eprev.get());
//                        auto insert2 = inserted_edges.emplace(i.enext.get());
//                        if (insert1.second) ret.add_value(i.eprev);
//                        if (insert2.second) ret.add_value(i.enext);
//                }
//        }
//        auto rtree = HMCont2D::Assembler::ETree(ret);
//        //match directions of all contours in rtree with directions in cont
//        //!!! This could change directions of closed contours of
//        //    ExtendedTree which are now directed according to its nested structure.
//        //    However, we hope that nested structure of ExtendedTree
//        //    will not be used later in the algorithm and this will not matter.
//        for (auto c: rtree.all_contours()){
//                assert(c->size() > 1);  //since we compose contours from point neighbors
//                Point* p1 = c->first();
//                Point* p2 = c->next_point(p1);
//                HMCont2D::Contour* parcont = cont.get_contour(p1);
//                Point* p3 = parcont->next_point(p1);
//                if (p3 != p2) c->Reverse();
//        }
//        return rtree;
//}

HMCont2D::Container<HMCont2D::ExtendedTree> do_bnd_partition(
		HMCont2D::ExtendedTree& segments, std::map<Point*, double>& bw,
		std::set<Point*>& needp){
	HMCont2D::ExtendedTree ret;
	HMCont2D::PCollection pstore;
	for (auto c: segments.all_contours()){
		std::map<double, double> bw1;
		for (auto it: bw){
			if (c->contains_point(it.first)){
				auto z = c->coord_at(*it.first);
				bw1[std::get<1>(z)] = it.second;
			}
		}
		vector<Point*> keep_p;
		for (auto p: c->all_points()){
			if (needp.find(p) != needp.end()) keep_p.push_back(p);
		}
		HMCont2D::Contour part1 = HMCont2D::Algos::WeightedPartition(
				bw1, *c, pstore, keep_p);
		ret.AddContour(part1);
	}
	return HMCont2D::Container<HMCont2D::ExtendedTree>::DeepCopy(ret);
}

void subpath_substitute(HMCont2D::ContourTree& cont, HMCont2D::ExtendedTree& repart){
	for (auto c: repart.all_contours()){
		Point* a1 = c->first();
		Point* a2 = c->last();
		Point* b1 = HMCont2D::ECollection::FindClosestNode(cont, *a1);
		Point* b2 = HMCont2D::ECollection::FindClosestNode(cont, *a2);
		assert(Point::dist(*a1, *b1) < geps && Point::dist(*a2, *b2) < geps);
		HMCont2D::Contour* subc = cont.get_contour(b1);
		assert(subc->contains_point(b1) && subc->contains_point(b2));
		auto ap = subc->ordered_points();
		int i1 = std::find(ap.begin(), ap.end(), b1) - ap.begin();
		int i2 = std::find(ap.begin(), ap.end(), b2) - ap.begin();
		if (subc->is_closed()){
			std::rotate(subc->data.begin(), subc->data.begin() + i1, subc->data.end());
			i2 -= i1;
			i1 -= i1;
			if (i2<=i1) i2+=subc->size();
		}
		//cut edges
		subc->data.erase(subc->data.begin() + i1, subc->data.begin() + i2);
		//paste edges
		subc->data.insert(subc->data.begin()+i1, c->data.begin(), c->data.end());
		//substitute equal points pointers
		for (auto e: *subc){
			if (e->pstart == a1) e->pstart = b1;
			if (e->pstart == a2) e->pstart = b2;
			if (e->pend == a1) e->pend = b1;
			if (e->pend == a2) e->pend = b2;
		}
	}
	cont.ReloadEdges();
}

void substitute_bw(HMCont2D::ExtendedTree& oldt, HMCont2D::ExtendedTree& newt,
		std::map<Point*, double>& bw){
	//1) remove unused points
	for (auto c: oldt.all_contours()){
		auto ap = c->ordered_points();
		//keep start/end points since they were not changed in subpath_substitute routine
		for (int i=1; i<(int)ap.size()-1; ++i){
			auto fnd = bw.find(ap[i]);
			if (fnd != bw.end()) bw.erase(fnd);
		}
	}
	//2) add partition points
	for (auto c: newt.all_contours()){
		auto ap = c->ordered_points();
		vector<double> lens = HMCont2D::ECollection::ELengths(*c);
		//first last points
		if (c->is_closed()) bw[ap[0]] = std::max(lens[0], lens.back());
		else{
			if (bw.find(ap[0]) != bw.end()) bw[ap[0]] = std::max(lens[0], bw[ap[0]]);
			else bw[ap[0]] = lens[0];
			if (bw.find(ap.back()) != bw.end()) bw[ap.back()] = std::max(lens.back(), bw[ap.back()]);
			else bw[ap.back()] = lens.back();
		}
		//internal
		for (int i=1; i<(int)ap.size()-1; ++i) bw[ap[i]] = std::max(lens[i-1], lens[i]);
	}
}

void substitute_bnded(HMCont2D::ExtendedTree& oldt, HMCont2D::ExtendedTree& newt,
		HMCont2D::ECollection& bnded){
	//1) remove old
	std::set<int> bade;
	for (auto c: oldt.all_contours()){
		for (auto e: c->data) bade.insert(bnded.get_index(e.get()));
	}
	bnded.RemoveAt(bade);
	//2) add new
	for (auto e: newt) bnded.add_value(e);
}

void simplify_bnd_edges(HMCont2D::ContourTree& cont, std::map<Point*, double>& bw,
                HMCont2D::ECollection& bnded, HMCont2D::PCollection& pstore, double angle0){
	if (ISZERO(angle0)) angle0 = 1e-6;
	else if (angle0<0) return;
	//needed points from bnded collection
	std::set<Point*> needed_points = find_needed_points(cont, bnded, angle0);
	//segments for partition
	HMCont2D::ExtendedTree segments = insignificant_segments(cont, needed_points);
	if (segments.cont_count() == 0) return;
	//do partition
	HMCont2D::Container<HMCont2D::ExtendedTree> repart = do_bnd_partition(segments, bw, needed_points);
	//substitute subpaths in cont
	subpath_substitute(cont, repart);
	pstore = repart.pdata;   //save all points since they are used by cont
	//substitute bnded
	substitute_bnded(segments, repart, bnded);
	//substitute bw
	substitute_bw(segments, repart, bw);
}

void simplify_source_edges(HMCont2D::ContourTree& cont, const HMCont2D::PCollection& srcpnt,
		std::map<const HMCont2D::Edge*, int>& edtypes){
	std::vector<const Point*> bad_points;
	for (auto c: cont.nodes){
		auto cinfo = c->ordered_info(); cinfo.resize(cinfo.size()-1);
		for (auto& ci: cinfo){
			if (edtypes[ci.eprev.get()] == 1 && edtypes[ci.enext.get()] == 1){
				auto dst = HMCont2D::PCollection::FindClosestNode(srcpnt, *ci.p);
				if (!ISZERO(sqrt(std::get<2>(dst)))) bad_points.push_back(ci.p);
			}
		}
	}
	cont.RemovePoints(bad_points);
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
> BufferGrid::boundary_info(bool preserve_bp, double angle0) const{
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
	simplify_source_edges(cont, scont.pdata, edtypes);

	//simplify bnd edges if needed:
	//  for connected sequence of bnd_edges changes last point of first edge
	//  and deletes all others from contour tree
	//if (!preserve_bp) simplify_bnd_edges(cont, edtypes, angle0);

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
	if (inner_edges.size() == 0 || outer_edges.size() == 0){
		//if no edges of certain type just weight all possible points weights by their distances
		std::map<Point*, double> allp;
		for (auto d: bw) if (d.second>0) allp.insert(d);
		//leave less then 100 nodes for performance reasons.
		//!!! Here a better approach would be to divide all allp points bounding box
		//into subsquares and get a point from each of them.
		if (allp.size() > 200){
			int k = allp.size() / 100;
			int i = 0;
			auto it = allp.begin();
			while (it!=allp.end()){
				if (i % k != 0) it = allp.erase(it);
				else ++it;
				++i;
			}
		}
		//go on weighting
		for (auto& d: bw) if (d.second<=0){
			d.second = 0;
			double dsum = 0;
			for (auto p: allp){
				double dist = Point::dist(*d.first, *p.first);
				dsum += 1.0/dist;
				d.second += p.second / dist;
			}
			d.second /= dsum;
		}
	} else
	for (auto& d: bw) if (d.second<0){
		//weight closest inner and closest outer point.
		//auto d1 = HMCont2D::ECollection::FindClosestEdge(inner_edges, *d.first);
		//auto d2 = HMCont2D::ECollection::FindClosestEdge(outer_edges, *d.first);
		//double dist1 = std::get<1>(d1), dist2 = std::get<1>(d2);
		//double len1 = std::get<0>(d1)->length(), len2 = std::get<0>(d2)->length();
		auto pc1 = HMCont2D::ECollection::FindClosestNode(inner_edges, *d.first);
		auto pc2 = HMCont2D::ECollection::FindClosestNode(outer_edges, *d.first);
		double dist1 = Point::dist(*pc1, *d.first), dist2 = Point::dist(*pc2, *d.first);
		double len1 = bw[pc1], len2 = bw[pc2];
		assert(len1>0 && len2>0);
		double w1 = dist1/(dist1 + dist2);
		d.second = (1 - w1)*len1 + w1*len2;
	}

	//simplify bnd edges if needed:
	//  for connected sequence of bnd_edges changes last point of first edge
	//  and deletes all others from contour tree
	if (!preserve_bp) simplify_bnd_edges(cont, bw, bnd_edges, apoints, angle0);

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
	//auto flt = cont_orig.filter_points_i(nocong);
	//if (std::get<1>(flt).size() != nocong.size())
		//throw std::runtime_error("Buffer grid does not match original");

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
