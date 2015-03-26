#include "intrusion.h"
#include "assert.h"

//find a cell containing point
int GridForIntrusion::find_cell(const Point* p){
	//it is supposed that this procedure will be called rarely
	//Hence implementing naive algorithm
	bool hint = true;
	for (int i=0; i<n_cells(); ++i){
		auto pc = cells[i]->get_contour();
		if (pc.is_inside(*p, &hint) != OUTSIDE) return i;
	}
	throw std::runtime_error("Failed to find point in any cell");
}

//delete a cell which contains collection of grids
//adds this collection supporting cells non-overlapping
//and single connection property
void GridForIntrusion::intrude(int icell, const vector<GridGeom*>& subs){
	ContoursCollection cc;
	cc.add_contour(get_cell(icell)->get_contour());
	for (auto& s: subs){
		auto c = s->get_contours();
		assert(c.size() == 1);
		cc.add_contour(c[0]);
	}
#ifndef NDEBUG
	assert(cc.get_level(0) == 0);
	for (int i=1; i<cc.n_cont(); ++i) assert(cc.get_level(i) == 1);
#endif
	//remove target cell
	aa::remove_entries(cells, {icell});
	//add all points and cells from subs
	std::map<const GridPoint*, const GridPoint*> oldnew_points;
	for (auto& s: subs){
		for (int i=0; i<s->n_points(); ++i){
			auto p = s->get_point(i);
			auto newp = aa::add_shared(points, GridPoint(*p));
			oldnew_points[p] = newp;
		}
		for (int i=0; i<s->n_cells(); ++i){
			auto newc = aa::add_shared(cells, Cell());
			auto oldc = s->get_cell(i);
			for (int j=0; j<oldc->dim(); ++j){
				auto newp = const_cast<GridPoint*>(oldnew_points[oldc->get_point(j)]);
				add_point_to_cell(newc, newp);
			}
		}
	}
	vector<vector<const Point*>> filling = noncross_filling(cc, 0);
	for (auto& f: filling){
		auto newc = aa::add_shared(cells, Cell());
		for (auto& p: f){
			auto gp = static_cast<const GridPoint*>(p);
			auto fnd = oldnew_points.find(gp);
			GridPoint* p1;
			if (fnd == oldnew_points.end()) p1 = const_cast<GridPoint*>(gp);
			else p1 = const_cast<GridPoint*>(fnd->second);
			add_point_to_cell(newc, p1);
		}
	}
}

//all childs should lie inside one of the parent cells
void GridForIntrusion::intrude_grids(GridGeom& parent, const vector<GridGeom*>& childs){
	std::map<int, vector<GridGeom*>> cell_childs;
	GridForIntrusion ig(parent);
	for (auto ch: childs){
		int icell = ig.find_cell(ch->get_point(0));
		if (cell_childs.find(icell) == cell_childs.end())
			cell_childs[icell] = vector<GridGeom*>();
		cell_childs[icell].push_back(ch);
	}
	//add cells from highest index to lowest no to break cell indexing
	for (auto it=cell_childs.rbegin(); it!=cell_childs.rend(); ++it){
		ig.intrude(it->first, it->second);
	}
	ig.swap_back(parent);
}

// ================= Contours collection filling
//noncross_filling procedures
namespace{

void noncross_impl(vector<ContoursCollection>& data){
	//----- recursion tail
	int in = -1;
	for (int i=0; i<data.size(); ++i)
		if (data[i].n_cont() > 1) { in = i; break; }
	if (in<0) return;

	//---- contour processing
	ContoursCollection& cc=data[in];
	PContour main=cc.contour(0), sec=cc.contour(1);
	
	//secondary contour points candidates
	auto keycomp = [](const Point* p1, const Point* p2){return *p1<*p2;};
	std::set<const Point*, decltype(keycomp)> src_cand(keycomp);
	for (int i=0; i<sec.n_points(); ++i) src_cand.insert(sec.get_point(i));

	//main contour points candidates
	std::set<const Point*> tar_cand;
	for (int i=0; i<main.n_points(); ++i) tar_cand.insert(main.get_point(i));
	
	//borders which are forbidden to cross
	std::vector<std::pair<const Point*, const Point*>> nocross;
	for (int i=0; i<cc.n_cont(); ++i){
		auto cont = cc.get_contour(i);
		for (int j=0; j<cont->n_points(); ++j){
			nocross.push_back(std::make_pair(cont->get_point(i), cont->get_point(i+1)));
		}
	}
	auto is_good_point = [&](const Point* p1, const Point* p2)->bool{
		double ksieta[2];
		for (auto& nc: nocross){
			SectCross(*nc.first, *nc.second, *p1, *p2, ksieta);
			if (ksieta[1]>geps && ksieta[1]<1-geps) return false;
		}
		return true;
	};
	
	//finding two connections
	//  c1.first = c1p0 - point on the internal contour
	//  c1.second = c1p1 - point on the outer contour
	//  c2.first = c2p0 - point on the internal contour
	//  c2.second = c2p1 - point on the outer contour
	std::pair<int, int> c1, c2;
	const Point *c1p0=0, *c1p1=0, *c2p0=0, *c2p1=0;
	
	auto try_c2p0 = [&]()->bool{
		//known c1p0, c1p1, c2p0
		std::map<double, const Point*> meas;
		for (auto& p: tar_cand) meas[Point::meas(*c2p0, *p)]=p;
		for (auto& v: meas){
			c2p1 = v.second;
			if (is_good_point(c2p0, c2p1)) break;
			else c2p1 = 0;
		}
		return c2p1 != 0;
	};

	auto try_c1p1 = [&]()->bool{
		//known c1p0, c1p1
		tar_cand.erase(c1p1);
		for (auto it = src_cand.rbegin(); it!=src_cand.rend(); ++it){
			if (*it == c1p0) continue;
			c2p0 = *it;
			if (try_c2p0()) break;
			else c2p0 = 0;
		}
		tar_cand.insert(c1p1);
		return c2p0!=0;
	};

	auto try_c1p0 = [&]()->bool{
		//fill others from known c1p0
		src_cand.erase(c1p0);  //remove point from candidates list
		std::map<double, const Point*> meas;
		for (auto& p: tar_cand) meas[Point::meas(*c1p0, *p)]=p;
		for (auto& v: meas){
			c1p1 = v.second;
			if (is_good_point(c1p0, c1p1)){
				nocross.push_back(std::make_pair(c1p0, c1p1));
				if (try_c1p1()) break;
				else{
					nocross.resize(nocross.size()-1);
					c1p1 = 0;
				}
			} else c1p1 = 0;
		}
		src_cand.insert(c1p0); //paste point back to candidates list
		return c1p1 != 0;
	};

	for (auto it=src_cand.begin(); it!=src_cand.end(); ++it){
		c1p0 = *it;
		if (try_c1p0()) break;
		else c1p0 = 0;
	}
	if (c1p0 == 0) throw std::runtime_error("unable to fill contours collections");

	auto get_index = [&](const Point* p, PContour* c)->int{
		for (int i=0; i<c->n_points(); ++i){
			if (c->get_point(i) == p) return i;
		}
		throw std::runtime_error("Unable to find point in a contour");
	};

	c1.first = get_index(c1p0, &sec); c1.second = get_index(c1p1, &main);
	c2.first = get_index(c2p0, &sec); c2.second = get_index(c2p1, &main);
	
	//assembling two main contours for next recursion step
	PContour main2, main3;
	//inner contour points
	bool flag = true;
	for (int i=0; i<sec.n_points()+1; ++i){
		auto p = const_cast<Point*>(sec.get_point(i+c1.first));
		(flag) ? main2.add_point(p) : main3.add_point(p);
		if (p == sec.get_point(c2.first)){
			main3.add_point(p);
			flag = false;
		}
	}
	//outer contour points
	flag = true;
	for (int i=0; i<main.n_points()+1; ++i){
		auto p = const_cast<Point*>(main.get_point(i+c2.second));
		(flag) ? main2.add_point(p) : main3.add_point(p);
		if (p == main.get_point(c1.second)){
			main3.add_point(p);
			flag = false;
		}
	}

	//assembling contours collections
	ContoursCollection cc2, cc3;
	cc2.add_contour(main2); cc3.add_contour(main3);
	bool hint = true;
	for (int i=2; i<cc.n_cont(); ++i){
		auto cs = cc.get_contour(i);
		if (main2.is_inside(*cs->get_point(0), &hint) == INSIDE) cc2.add_contour(*cs);
		else cc3.add_contour(*cs);
	}

	//---- next call
	data[in] = cc2; data.push_back(cc3);
	noncross_impl(data);
}

}//namespace

vector<vector<const Point*>> GridForIntrusion::noncross_filling(const ContoursCollection& col, int i){
	vector<vector<const Point*>> ret;
	ContoursCollection subcol;
	subcol.add_contour(*col.get_contour(i));
	for (auto c: col.get_childs(i)) subcol.add_contour(*c);
	std::vector<ContoursCollection> dt; dt.reserve(col.n_cont());
	dt.push_back(col);
	noncross_impl(dt);
	for (auto cc: dt){
		ret.push_back(vector<const Point*>());
		auto& v = ret.back();
		auto c=cc.get_contour(i);
		for (int i=0; i<c->n_points(); ++i) v.push_back(c->get_point(i));
	}
	return ret;
}

