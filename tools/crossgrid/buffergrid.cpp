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

	//2) check if cont and main intersection meets requirements for building the buffer
	{
		//if any point the main grid lies strongly within the contour -> proceed with the procedure
		auto fnd1 = std::find_if(meas.begin(), meas.end(), [](double x){ return x>geps2; });
		if (fnd1 != meas.end()) goto proceed_proc;
	}
	{
		//if less then two points of main lie on the cont -> stop procedure
		int bndcount = std::count_if(meas.begin(), meas.end(), [](double x){ return fabs(x)<geps2; });
		if (bndcount<2) goto stop_proc;
	}
	{
		//if any boundary main edge lies on cont -> proceed
		for (auto e: main.get_edges()) if (e.is_boundary()){
			if (fabs(meas[e.p1])<geps2 && fabs(meas[e.p2])<geps2) goto proceed_proc;
		}
	}
	//stop/proceed points
	stop_proc: return;
	proceed_proc:;

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

//boundary_info helper functions
namespace{

std::pair<const Point*, double> closest_point(const Point* p, std::set<const Point*>& src){
	const Point* resp = NULL; double resd = 1e200;
	for (auto sp: src){
		double meas = Point::meas(*p, *sp);
		if (meas < resd){
			resd = meas;
			resp = sp;
		}
	}
	return std::make_pair(resp, sqrt(resd));
};

void clean_noncorner_points(PContour& c, std::set<const Point*>& src, const PContour* cont = NULL){
	std::set<int> bad_pts;
	for (int i=0; i<c.n_points(); ++i){
		const Point* p = c.get_point(i);
		//if corner point -> ignore
		if (c.is_corner_point(i)) continue;

		//if coinsides with cont point -> ignore
		if (cont != NULL){
			for (int j=0; j<cont->n_points(); ++j){
				if (*p == *cont->get_point(j)){
					p = NULL;
					break;
				}
			}
			if (p == NULL) continue;
		}


		//if lies in src -> delete from c and src
		auto fnd = src.find(p);
		if (fnd != src.end()){
			src.erase(fnd);
			bad_pts.insert(i);
		}
	}
	c.delete_by_index(bad_pts);
}

void boundary_points_dist(const PContour& c, std::map<const Point*, double>& d, std::set<const Point*>& bp){
	//types: 1 for bp points, 0 for others
	vector<int> types(c.n_points(), 0);
	for (int i=0; i<c.n_points(); ++i){
		if (bp.find(c.get_point(i)) != bp.end()) types[i] = 1;
	}
	
	//loop over contour points
	for (int i=0; i<c.n_points(); ++i){
		const Point* p = c.get_point(i);
		int inext = (i==c.n_points()-1) ? 0 : i+1;
		int iprev = (i==0) ? c.n_points()-1 : i-1;

		if (types[i] == 1) d[p] = -1;
		else {
			if (types[iprev] == 1 && types[inext] == 1){
				d[p] = -1;
			} else if (types[iprev] == 1 && types[inext] == 0){
				d[p] = Point::dist(*p, *c.get_point(inext));
			} else if (types[inext] == 1 && types[iprev] == 0){
				d[p] = Point::dist(*p, *c.get_point(iprev));
			}
		}
	}
}

void inner_outer_points_dist(const PContour& c, std::map<const Point*, double>& d,
		std::set<const Point*>& bp1, std::set<const Point*>& bp2){
	//types: 1 for bp1 points, 2 for bp2
	vector<int> types(c.n_points(), 0);
	for (int i=0; i<c.n_points(); ++i){
		if (bp1.find(c.get_point(i)) != bp1.end()) types[i] = 1;
		else if (bp2.find(c.get_point(i)) != bp2.end()) types[i] = 2;
	}
	
	//loop over contour points
	for (int i=0; i<c.n_points(); ++i){
		const Point* p = c.get_point(i);
		int inext = (i==c.n_points()-1) ? 0 : i+1;
		int iprev = (i==0) ? c.n_points()-1 : i-1;

		int tpi = types[i], tpip = types[inext], tpim = types[iprev];

		if (tpi < 1 || tpip < 1 || tpim < 1) continue; //this was treated in boundary_points_dist
		else if (tpi == tpip == tpim) continue; //normal situation
		else if (tpi == tpip && tpi != tpim){
			d[p] = Point::dist(*p, *c.get_point(inext));
		}
		else if (tpi != tpim && tpi == tpim){
			d[p] = Point::dist(*p, *c.get_point(iprev));
		}
	}
}

}//namespace


std::tuple<
	ContoursCollection,
	std::vector<double>
> BufferGrid::boundary_info(bool preserve_true_bp) const{
	//return value
	std::tuple<ContoursCollection, vector<double>> ret;
	auto& cc = std::get<0>(ret);
	auto& lc = std::get<1>(ret);

	//nodes origin info
	auto origin = build_bedges();
	auto& inner_bp = std::get<0>(origin); //points on the source contour
	auto& outer_bp = std::get<1>(origin); //grid points outside bufferzone
	auto& true_bp = std::get<2>(origin);  //boundary grid points

	//if there is no outer and boundary points: return only source contour
	if (outer_bp.size() == 0 && true_bp.size() == 0){
		cc.add_contour(source_cont);
		lc = cc.get_contour(0)->chdist();
		return ret;
	}

	// --- build real contour
	//1) initial contour
	std::vector<PContour> real_contour = get_contours();

	//2) delete points which lie on the segments of the source contour and
	//are not its edge points. Remove points from c and from inner_bp.
	for (auto& c: real_contour) clean_noncorner_points(c, inner_bp, &source_cont);

	//3) assemble distances into map
	std::map<const Point*, double> dist;
	for (auto& c: real_contour){
		//calculatate distances for contour points
		auto cdist = c.chdist();
		for (int i=0; i<c.n_points(); ++i)
			dist[c.get_point(i)] = cdist[i];
		//treat distances for true boundary points. it should:
		//	- equal -1 by itself
		//	- not affect adjecent points distances
		boundary_points_dist(c, dist, true_bp);
		//treat situations when left side of point is inner and right side -- outer.
		//	- equal its distance to -1
		inner_outer_points_dist(c, dist, inner_bp, outer_bp);
	}
	
	//4) delete all non corner boundary points if necessary
	if (!preserve_true_bp){
		for (auto& c: real_contour) clean_noncorner_points(c, true_bp);
	}


	//5) calculate distances for points with illegal distances (== -1).
	//as the linear weighted combination of closest inner and closest outer points
	//distances
	auto weight = [&](const Point* p)->double{
		auto d1 = closest_point(p, inner_bp);
		auto d2 = closest_point(p, outer_bp);
		if (d1.first == 0) return dist[d2.first];
		if (d2.first == 0) return dist[d1.first];
		double w1 = d1.second/(d1.second + d2.second);
		return (1 - w1)*dist[d1.first] + w1*dist[d2.first];
	};
	for (auto& v: dist) if (v.second < 0){
		v.second = weight(v.first);
	}

	//6) assemble result
	for (auto& c: real_contour){
		//add contour and distances to result
		cc.add_contour(c);
	}
	for (int i=0; i<cc.n_cont(); ++i){
		auto cont = cc.get_contour(i);
		for (int i=0; i<cont->n_points(); ++i)
			lc.push_back(dist[cont->get_point(i)]);
	}

	return ret;
}

std::tuple<
	ContoursCollection,
	std::vector<double>
> BufferGrid::boundary_info2(bool preserve_true_bp) const{
	//return value
	std::tuple<ContoursCollection, vector<double>> ret;
	auto& cc = std::get<0>(ret);
	auto& lc = std::get<1>(ret);

	//nodes origin info
	auto origin = build_bedges();
	auto& inner_bp = std::get<0>(origin); //points on the source contour
	auto& outer_bp = std::get<1>(origin); //points grid points outside bufferzone
	auto& true_bp = std::get<2>(origin);  //boundary grid points

	//if there is no outer points: set the furthest true_bp point as outer
	//if there is no outer and boundary points: return only source contour
	std::shared_ptr<Point> dummy_outer;
	if (outer_bp.size()==0){
		if (true_bp.size()>0){
			double maxmeas = 0;
			const Point* furthp = 0;
			for (auto& ptr: true_bp){
				double minmeas = 1e100;
				for (auto& pin: inner_bp){
					double m = Point::meas(*ptr, *pin);
					if (m<minmeas) minmeas = m;
				}
				if (minmeas>maxmeas){
					maxmeas = minmeas;
					furthp = ptr;
				}
			}
			//dummy_outer.reset(new Point(furthp->x, furthp->y));
			//outer_bp.insert(dummy_outer.get());
			true_bp.erase(true_bp.find(furthp));
			outer_bp.insert(furthp);
		} else {
			cc.add_contour(source_cont);
			lc = source_cont.chdist();
			return ret;
		}
	}

	// --- build real contour
	std::vector<PContour> real_contour = get_contours();
	//delete points which lie on the sections of the source contour
	//delete all boundary points if necessary
	for (auto& c: real_contour){
		std::set<int> bad_points;
		for (int i=0; i<c.n_points(); ++i){
			if (!c.is_corner_point(i)){
				auto p = c.get_point(i);
				if (inner_bp.find(p)!=inner_bp.end()){
					//find if p is the edge of source contour
					for (int j=0; j<source_cont.n_points(); ++j){
						if (*p == *source_cont.get_point(j)){
							p=0; break;
						}
					}
					if (p!=0) bad_points.insert(i);
				} else if (!preserve_true_bp && true_bp.find(p)!=true_bp.end()){
					//delete all boundary points
					bad_points.insert(i);
				}
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
	//filter points origin to preserve only points from dist_dict
	auto filter_origin = [&dist_dict](std::set<const Point*>& st){
		auto it = st.begin();
		while (it!=st.end()){
			if (dist_dict.find(*it)==dist_dict.end()){
				auto it2 = it++;
				st.erase(it2);
			} else ++it;
		}
	};
	filter_origin(inner_bp);
	filter_origin(outer_bp);
	filter_origin(true_bp);
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
		double f1 = *dist_dict[closest_inner.first];
		double f2 = *dist_dict[closest_outer.first];
		double d1 = closest_inner.second, d2 = closest_outer.second;
		*dist_dict[p] = (f1*d2+f2*d1)/(d1+d2);
	}

	//assemble contours collection and return
	for (auto& c: real_contour) cc.add_contour(c);
	return ret;
}

void BufferGrid::update_original() const{
	//1) delete all original cells
	std::set<int> _ds;
	aa::Cfill_container(orig_cells, _ds, [](const Cell* p){ return p->get_ind(); });
	aa::remove_entries(orig->cells, _ds);

	//2) backup original boundary points
	std::set<const GridPoint*> bpts_orig = orig->get_bnd_points();

	//3) add all buffer cells and nodes
	//to the end of original data arrays
	orig->add_data(*this);

	//4) get all new boundary points
	std::set<const GridPoint*> bpts_tmp = orig->get_bnd_points();
	std::set<const GridPoint*> bpts_new;
	for (auto x: bpts_tmp)
		if (bpts_orig.find(x)==bpts_orig.end())
			bpts_new.insert(x);

	//5) find all congruent contour points
	//original point index -> buffer point index
	std::map<int, int> cong_p;
	for (auto i=bpts_orig.begin(); i!=bpts_orig.end(); ++i){
		for (auto j=bpts_new.begin(); j!=bpts_new.end(); ++j){
			if (**i == **j){
				cong_p.emplace((*i)->get_ind(), (*j)->get_ind());
				break;
			}
		}
	}

	//6) remove points of the original grid which lie on the source contour.
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

	//7) merge congruent points and remove bad points
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

	//6) delete unused
	orig->delete_unused_points();
	orig->set_indicies();
}

