#include <map>
#include "addalgo.hpp"
#include "buffergrid.h"
#include "procgrid.h"
#include "wireframegrid.h"
#include "cont_assembler.hpp"
#include "constructor.hpp"
#include "algos.hpp"
#include "contclipping.hpp"
#include "cont_partition.hpp"
namespace{

//this is a temporary function which will be deleted after
//crossgrid::Contour class will be completely substituted by HM2D classes
HM2D::EdgeData contour_to_hm(const Contour& c){
	std::vector<Point> pnt;
	for (int i=0; i<c.n_points(); ++i) pnt.push_back(*c.get_point(i));
	return HM2D::Contour::Constructor::FromPoints(pnt, true);
}

HM2D::EdgeData get_contact(HM2D::Contour::Tree& gridcont, HM2D::EdgeData& source){
	PtsGraph g(gridcont.alledges()), s(source);
	PtsGraph over = PtsGraph::overlay(g, s);
	//leave only those lines which centers lie on source
	//and inside gridcont
	vector<Point> cp = over.center_line_points();
	vector<int> srt = HM2D::Contour::Algos::SortOutPoints(gridcont, cp);
	std::set<int> bad_lines;
	for (int i=0; i<cp.size(); ++i) if (srt[i] != INSIDE) bad_lines.insert(i);
	over.exclude_lines(bad_lines);
	return over.toecollection();
}

//calculates measures to contact lines for points in pt.
//if points lies without bounding box of given contours + buf sets 2*buf*buf as its measure
vector<double> build_measures(HM2D::EdgeData& src, const ShpVector<GridPoint>& pt){
	vector<double> ret;
	for (auto p: pt){
		auto r = HM2D::FindClosestEdge(src, *p);
		double& av = std::get<1>(r);
		ret.push_back(av*av);
	}
	return ret;
}

//for points wich lie outside given area multiplies its measure by -1
void negative_to_outside(vector<double>& meas, const HM2D::EdgeData& area, const ShpVector<GridPoint>& pt){
	vector<Point> p2;
	for (auto p: pt) p2.push_back(*p);
	vector<int> srt = HM2D::Contour::Algos::SortOutPoints(area, p2);
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
	HM2D::Contour::Tree cross;
	if (source_cont.area() > 0) cross = HM2D::Contour::Clip::Difference(treeorig, hmcont);
	else cross = HM2D::Contour::Clip::Intersection(treeorig, hmcont);
	HM2D::Contour::Clip::Heal(cross);
	double domarea = cross.area();
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
		int pnts_with_zero_m = 0;
		for (int j=0; j<c->dim(); ++j){
			double& m = meas[c->get_point(j)->get_ind()];
			//point on contour are ambiguous
			if (fabs(m)<geps2){
				++pnts_with_zero_m;
				if (j!=c->dim()-1){
					is_good = true; continue;
				} else if (pnts_with_zero_m < c->dim()){
					is_good =  true;
				}else{
					//all points lie on source contour. we have to check intersection
					auto ccont = GGeom::Info::CellContour(main, c->get_ind());
					auto intersect = HM2D::Contour::Clip::Union(cross, ccont);
					HM2D::Contour::Clip::Heal(intersect);
					double unarea = intersect.area();
					if (fabs(unarea-domarea)>100*geps*geps) is_good = false;
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

//void simplify_bnd_edges(HM2D::Contour::Tree& cont,
//                std::map<const HM2D::Edge*, int>& edtypes, double angle0){
//        if (ISZERO(angle0)) angle0 = 1e-6;
//        //assemble needed points because there could be cases when a point
//        //is needed from one contour of 'cont' and not needed from the other
//        std::set<const Point*> needed_points;
//        for (auto c: cont.nodes){
//                vector<HM2D::EdgeData::PInfo> info = c->ordered_info();
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

std::set<Point*> find_needed_points(HM2D::Contour::Tree& cont, HM2D::EdgeData& bnded, double a0){
	a0 *= (M_PI/180);
	std::set<Point*> ret;
	for (auto c: cont.nodes){
		auto info = HM2D::Contour::OrderedInfo(c->contour);
		for (auto i: info){
			if (!i.eprev || !i.enext ||
			    !HM2D::Contains(bnded, i.eprev.get()) ||
			    !HM2D::Contains(bnded, i.enext.get()) ||
			    fabs(Angle(*i.pprev, *i.p, *i.pnext)-M_PI)>a0){
				ret.insert(i.p.get());
			}
		}
	}
	return ret;
}

HM2D::Contour::Tree insignificant_segments(HM2D::Contour::Tree& cont, std::set<Point*>& needp){
	std::set<HM2D::Edge*> inserted_edges;
	HM2D::EdgeData ret;
	for (auto c: cont.nodes){
		auto info = HM2D::Contour::OrderedInfo(c->contour);
		for (auto i: info) if (needp.find(i.p.get()) == needp.end()){
			auto insert1 = inserted_edges.emplace(i.eprev.get());
			auto insert2 = inserted_edges.emplace(i.enext.get());
			if (insert1.second) ret.push_back(i.eprev);
			if (insert2.second) ret.push_back(i.enext);
		}
	}
	auto rtree = HM2D::Contour::Tree::Assemble(ret);
	//match directions of all contours in rtree with directions in cont
	for (auto c: rtree.nodes){
		assert(c->contour.size() > 1);  //since we compose contours from point neighbors
		Point* p1 = HM2D::Contour::First(c->contour).get();
		Point* p2 = c->contour[0]->sibling(p1).get();
		HM2D::EdgeData* parcont = &cont.find_node(p1)->contour;
		auto i3 = HM2D::Contour::PInfo(*parcont, p1);
		if (i3.pnext.get() != p2) HM2D::Contour::Reverse(c->contour);
	}
	return rtree;
}
//HM2D::Contour::Tree insignificant_segments(HM2D::Contour::Tree& cont, HM2D::EdgeData& bnded, double a0){
//        std::set<HM2D::Edge*> inserted_edges;
//        HM2D::EdgeData ret;
//        for (auto c: cont.nodes){
//                vector<HM2D::EdgeData::PInfo> info = c->ordered_info();
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
//        auto rtree = HM2D::Assembler::ETree(ret);
//        //match directions of all contours in rtree with directions in cont
//        //!!! This could change directions of closed contours of
//        //    ExtendedTree which are now directed according to its nested structure.
//        //    However, we hope that nested structure of ExtendedTree
//        //    will not be used later in the algorithm and this will not matter.
//        for (auto c: rtree.all_contours()){
//                assert(c->size() > 1);  //since we compose contours from point neighbors
//                Point* p1 = c->first();
//                Point* p2 = c->next_point(p1);
//                HM2D::EdgeData* parcont = cont.get_contour(p1);
//                Point* p3 = parcont->next_point(p1);
//                if (p3 != p2) c->Reverse();
//        }
//        return rtree;
//}

HM2D::Contour::Tree do_bnd_partition(
		HM2D::Contour::Tree& segments, std::map<Point*, double>& bw,
		std::set<Point*>& needp){
	HM2D::Contour::Tree ret;
	for (auto c: segments.nodes){
		std::map<double, double> bw1;
		for (auto it: bw){
			if (Contains(c->contour, it.first)){
				auto z = HM2D::Contour::CoordAt(c->contour, *it.first);
				bw1[std::get<1>(z)] = it.second;
			}
		}
		HM2D::VertexData keep_p;
		for (auto p: AllVertices(c->contour)){
			if (needp.find(p.get()) != needp.end()) keep_p.push_back(p);
		}
		HM2D::EdgeData part1 = HM2D::Contour::Algos::WeightedPartition(
				bw1, c->contour, keep_p);
		ret.AddContour(part1);
	}
	return HM2D::Contour::Tree::DeepCopy(ret);
}

void subpath_substitute(HM2D::Contour::Tree& cont, HM2D::Contour::Tree& repart){
	for (auto c: repart.nodes){
		Point* a1 = HM2D::Contour::First(c->contour).get();
		Point* a2 = HM2D::Contour::Last(c->contour).get();
		auto av = AllVertices(cont.alledges());
		auto tb1 = HM2D::FindClosestNode(av, *a1);
		auto tb2 = HM2D::FindClosestNode(av, *a2);
		auto b1 = av[std::get<0>(tb1)];
		auto b2 = av[std::get<0>(tb2)];
		assert(Point::dist(*a1, *b1) < geps && Point::dist(*a2, *b2) < geps);
		HM2D::EdgeData* subc = &cont.find_node(b1.get())->contour;
		assert(HM2D::Contains(*subc, b1.get()) && HM2D::Contains(*subc, b2.get()));
		auto ap = HM2D::Contour::OrderedPoints(*subc);
		int i1 = std::find(ap.begin(), ap.end(), b1) - ap.begin();
		int i2 = std::find(ap.begin(), ap.end(), b2) - ap.begin();
		if (HM2D::Contour::IsClosed(*subc)){
			std::rotate(subc->begin(), subc->begin() + i1, subc->end());
			i2 -= i1;
			i1 -= i1;
			if (i2<=i1) i2+=subc->size();
		}
		//cut edges
		subc->erase(subc->begin() + i1, subc->begin() + i2);
		//paste edges
		subc->insert(subc->begin()+i1, c->contour.begin(), c->contour.end());
		//substitute equal points pointers
		for (auto e: *subc){
			if (e->first().get() == a1) e->vertices[0] = b1;
			if (e->first().get() == a2) e->vertices[0] = b2;
			if (e->last().get() == a1) e->vertices[1] = b1;
			if (e->last().get() == a2) e->vertices[1] = b2;
		}
	}
}

void substitute_bw(HM2D::Contour::Tree& oldt, HM2D::Contour::Tree& newt,
		std::map<Point*, double>& bw){
	//1) remove unused points
	for (auto c: oldt.nodes){
		auto ap = HM2D::Contour::OrderedPoints(c->contour);
		//keep start/end points since they were not changed in subpath_substitute routine
		for (int i=1; i<(int)ap.size()-1; ++i){
			auto fnd = bw.find(ap[i].get());
			if (fnd != bw.end()) bw.erase(fnd);
		}
	}
	//2) add partition points
	for (auto c: newt.nodes){
		auto ap = HM2D::Contour::OrderedPoints(c->contour);
		vector<double> lens = HM2D::ELengths(c->contour);
		//first last points
		if (HM2D::Contour::IsClosed(c->contour)) bw[ap[0].get()] = std::max(lens[0], lens.back());
		else{
			if (bw.find(ap[0].get()) != bw.end()) bw[ap[0].get()] = std::max(lens[0], bw[ap[0].get()]);
			else bw[ap[0].get()] = lens[0];
			if (bw.find(ap.back().get()) != bw.end()) bw[ap.back().get()] = std::max(lens.back(), bw[ap.back().get()]);
			else bw[ap.back().get()] = lens.back();
		}
		//internal
		for (int i=1; i<(int)ap.size()-1; ++i) bw[ap[i].get()] = std::max(lens[i-1], lens[i]);
	}
}

void substitute_bnded(HM2D::Contour::Tree& oldt, HM2D::Contour::Tree& newt,
		HM2D::EdgeData& bnded){
	//1) remove old
	for (auto c: oldt.nodes)
	for (auto e: c->contour){
		auto fnd = std::find(bnded.begin(), bnded.end(), e);
		if (fnd != bnded.end()) bnded.erase(fnd);
	}
	//2) add new
	for (auto e: newt.alledges()) bnded.push_back(e);
}

void simplify_bnd_edges(HM2D::Contour::Tree& cont, std::map<Point*, double>& bw,
                HM2D::EdgeData& bnded, double angle0){
	if (ISZERO(angle0)) angle0 = 1e-6;
	else if (angle0<0) return;
	//needed points from bnded collection
	std::set<Point*> needed_points = find_needed_points(cont, bnded, angle0);
	//segments for partition
	HM2D::Contour::Tree segments = insignificant_segments(cont, needed_points);
	if (segments.nodes.size() == 0) return;
	//do partition
	HM2D::Contour::Tree repart = do_bnd_partition(segments, bw, needed_points);
	//substitute subpaths in cont
	subpath_substitute(cont, repart);
	//substitute bnded
	substitute_bnded(segments, repart, bnded);
	//substitute bw
	substitute_bw(segments, repart, bw);
}

void simplify_source_edges(HM2D::Contour::Tree& cont, const HM2D::VertexData& srcpnt,
		std::map<const HM2D::Edge*, int>& edtypes){
	for (auto c: cont.nodes){
		if (c->isopen()) continue;
		std::vector<const Point*> bad_points;
		auto cinfo = HM2D::Contour::OrderedInfo(c->contour);
		cinfo.resize(cinfo.size()-1);
		for (auto& ci: cinfo){
			if (edtypes[ci.eprev.get()] == 1 && edtypes[ci.enext.get()] == 1){
				auto ps = HM2D::FindClosestNode(srcpnt, *ci.p);
				if (*ci.p != *srcpnt[std::get<0>(ps)]) bad_points.push_back(ci.p.get());
			}
		}
		auto op = HM2D::Contour::OrderedPoints(c->contour);
		HM2D::VertexData op2;
		for (int i=0; i<op.size()-1; ++i)
		if (std::find(bad_points.begin(), bad_points.end(), op[i].get())
				== bad_points.end()){
			op2.push_back(op[i]);
		}
		if (op2.size() != op.size()-1){
			c->contour = HM2D::Contour::Assembler::Contour1(op2, true);
		}
	}
}

//extracts edges of 'from' which lie on 'where'
HM2D::EdgeData extract_edges(const HM2D::EdgeData& from,
		const HM2D::EdgeData& where){
	std::map<const Point*, double> pdist;
	for(auto p: AllVertices(from)){
		pdist[p.get()] = std::get<1>(HM2D::FindClosestEdge(where, *p));
	}
	HM2D::EdgeData ret;
	for(auto e: from){
		double dist1 = pdist[e->first().get()];
		double dist2 = pdist[e->last().get()];
		if (ISZERO(dist1) && ISZERO(dist2)) ret.push_back(e);
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
	HM2D::Contour::Tree,
	std::map<Point*, double>
> BufferGrid::boundary_info(bool preserve_bp, double angle0) const{
	//prepare return
	std::tuple<
		HM2D::Contour::Tree,
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
	HM2D::EdgeData inner_edges = extract_edges(cont.alledges(), scont);
	HM2D::EdgeData bnd_edges = extract_edges(cont.alledges(), grid_cont.alledges());
	HM2D::EdgeData outer_edges;
	std::map<const HM2D::Edge*, int> edtypes;
	for (auto& e: cont.alledges()) edtypes[e.get()] = 3;
	for (auto& e: inner_edges) edtypes[e.get()] = 1;
	for (auto& e: bnd_edges) edtypes[e.get()] = 2;

	//inner, bnd, outer edges set. Edges are shared with cont
	for (auto& ed: cont.alledges()){
		if (edtypes[ed.get()] == 3) outer_edges.push_back(ed);
	}
	assert(inner_edges.size()!=0 || outer_edges.size()!=0);

	// ================ cont purge procedures
	simplify_source_edges(cont, AllVertices(scont), edtypes);

	//simplify bnd edges if needed:
	//  for connected sequence of bnd_edges changes last point of first edge
	//  and deletes all others from contour tree
	//if (!preserve_bp) simplify_bnd_edges(cont, edtypes, angle0);

	//Edges deleted from cont still present in bnd_edges collection but have NULL data.
	//clear *_edges ecollections from NULL data edges for further weight calculations
	auto rmfun = [](shared_ptr<HM2D::Edge> e){return e->first() == nullptr;};
	inner_edges.resize(std::remove_if(inner_edges.begin(), inner_edges.end(), rmfun)-inner_edges.begin());
	outer_edges.resize(std::remove_if(outer_edges.begin(), outer_edges.end(), rmfun)-outer_edges.begin());
	bnd_edges.resize(std::remove_if(bnd_edges.begin(), bnd_edges.end(), rmfun)-bnd_edges.begin());

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
		auto oinfo = HM2D::Contour::OrderedInfo(c->contour);
		auto olens = HM2D::ELengths(c->contour);
		for (int i=0; i<c->contour.size(); ++i){
			double iprev = (i==0)?c->contour.size()-1:i-1;
			const HM2D::Edge *eprev=oinfo[i].eprev.get();
			const HM2D::Edge *enext=oinfo[i].enext.get();
			bw[oinfo[i].p.get()] = calc_weight(olens[iprev], olens[i], edtypes[eprev], edtypes[enext]);
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
		//auto d1 = HM2D::EdgeData::FindClosestEdge(inner_edges, *d.first);
		//auto d2 = HM2D::EdgeData::FindClosestEdge(outer_edges, *d.first);
		//double dist1 = std::get<1>(d1), dist2 = std::get<1>(d2);
		//double len1 = std::get<0>(d1)->length(), len2 = std::get<0>(d2)->length();
		auto iev = AllVertices(inner_edges);
		auto oev = AllVertices(outer_edges);
		auto rpc1 = HM2D::FindClosestNode(iev, *d.first);
		auto rpc2 = HM2D::FindClosestNode(oev, *d.first);
		auto pc1 = iev[std::get<0>(rpc1)].get();
		auto pc2 = oev[std::get<0>(rpc2)].get();
		double dist1 = Point::dist(*pc1, *d.first), dist2 = Point::dist(*pc2, *d.first);
		double len1 = bw[pc1], len2 = bw[pc2];
		assert(len1>0 && len2>0);
		double w1 = dist1/(dist1 + dist2);
		d.second = (1 - w1)*len1 + w1*len2;
	}
	//simplify bnd edges if needed:
	//  for connected sequence of bnd_edges changes last point of first edge
	//  and deletes all others from contour tree
	if (!preserve_bp) simplify_bnd_edges(cont, bw, bnd_edges, angle0);

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
