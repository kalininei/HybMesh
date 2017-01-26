#include "healgrid.hpp"
#include "algos.hpp"
#include "modgrid.hpp"
#include "treverter2d.hpp"
#include "cont_assembler.hpp"
#include "contabs2d.hpp"
#include "finder2d.hpp"
#include "contclipping.hpp"

using namespace HM2D;
using namespace HM2D::Grid;

namespace {

struct _ShiftSnapPreCalc{
	GridData* g;
	//HM2D::Contour::Tree gridbnd;
	std::map<double, shared_ptr<HM2D::Vertex>> contw;
	HM2D::VertexData cp;
	std::map<const Vertex*, double> bndw;
	HM2D::EdgeData bndeds;

	_ShiftSnapPreCalc(GridData& grid, const HM2D::EdgeData& cont,
			const VertexData& snap_nodes): g(&grid){
		auto cv = HM2D::AllVertices(cont);
		//snapping nodes
		for (auto p: snap_nodes){
			//try to search amoung vertices
			auto tfpnt = HM2D::Finder::ClosestPoint(cv, *p);
			Point* fpnt = cv[std::get<0>(tfpnt)].get();
			if (*fpnt == *p){
				p->set(*fpnt);
				continue;
			}
			//snap to edge
			auto fed = HM2D::Finder::ClosestEdge(cont, *p);
			HM2D::Edge* e = cont[std::get<0>(fed)].get();
			double w = std::get<2>(fed);
			p->set(Point::Weigh(*e->first(), *e->last(), w));
		}
		//gridbnd = GGeom::Info::Contour(grid);
		//all contour significant points weights
		cp = HM2D::Contour::CornerPoints(cont);
		for (auto p: cp){
			auto coord = HM2D::Contour::CoordAt(cont, *p);
			contw[std::get<1>(coord)] = p;
		}
		//copy one more time with w+1 for closed contours
		if (HM2D::Contour::IsClosed(cont)){
			auto it = contw.rbegin();
			while (it!=contw.rend()){
				contw[it->first + 1.0] = it->second;
				++it;
			}
		}
		bndeds = ECol::Assembler::GridBoundary(grid);
		//grid boundary points which lie on cont weights
		for (auto p: AllVertices(bndeds)){
			auto coord = HM2D::Contour::CoordAt(cont, *p);
			if (std::get<4>(coord)<geps) bndw[p.get()] = std::get<1>(coord);
		}
	}
	std::list<Point*> get_pts_between(double w1, double w2){
		std::list<Point*> ret;
		if (ISZERO(w2-w1)) return ret;
		if (w2<w1) w2+=1.0;
		assert(w1 < w2 && w2<2.0);
		for (auto it=contw.lower_bound(w1); it!=contw.end(); ++it){
			if (it->first>w1+geps && it->first<w2-geps){
				ret.push_back(it->second.get());
			}
			if (it->first>w2-geps) break;
		}
		return ret;
	}

	std::vector<
		std::tuple<Edge*, std::list<Point*>>
	>
	point_pairs(){
		std::vector<std::tuple<Edge*, std::list<Point*>>> ret;
		for (auto e: bndeds){
			Vertex *p1 = e->first().get(), *p2 = e->last().get();
			auto fnd1 = bndw.find(p1), fnd2 = bndw.find(p2);
			//if edge is not on contour ignore it
			if (fnd1==bndw.end() || fnd2==bndw.end()) continue;
			double w1 = fnd1->second, w2 = fnd2->second;
			std::list<Point*> ap = get_pts_between(w1, w2);
			if (ap.size() == 0) continue;
			Cell* cell = e->left.lock().get();
			ret.push_back(std::make_tuple(e.get(), ap));
		}
		return ret;
	}
};

}

void Algos::SnapToContour(GridData& grid, const HM2D::EdgeData& cont,
		const VertexData& snap_nodes){
	assert(Contour::IsContour(cont));
	//LeftCells reverter has index based algorithm. All new edges will be
	//added to the end of grid.vedges so it is save to used grid.vedges as
	//input data to the reverter.
	Contour::R::LeftCells grev(grid.vedges);
	Contour::R::ReallyDirect contrev(cont);

	auto proc = _ShiftSnapPreCalc(grid, cont, snap_nodes);
	aa::enumerate_ids_pvec(grid.vedges);
	for (auto& p: proc.point_pairs()){
		Edge* e = std::get<0>(p);
		vector<Point> padds;
		for (auto ap: std::get<1>(p)) padds.push_back(*ap);
		SplitEdge(grid, e->id, padds);
	}
}

namespace{
struct ShiftS{
	Vertex* from;
	Point* to;
	double dist2;
};
bool operator<(const ShiftS& a, const ShiftS& b){ return a.dist2<b.dist2;}
}

void Algos::ShiftToContour(GridData& grid, const HM2D::EdgeData& cont,
		const VertexData& snap_nodes){
	assert(Contour::IsContour(cont));
	Contour::R::LeftCells grev(grid.vedges);
	Contour::R::Clockwise contrev(cont, false);

	auto proc = _ShiftSnapPreCalc(grid, cont, snap_nodes);
	//get all possible shifts
	std::list<ShiftS> allshifts;
	for (auto p: proc.point_pairs()){
		Edge* e = std::get<0>(p);
		Vertex* p1 = e->first().get(), *p2 = e->last().get();
		auto& ap = std::get<1>(p);
		ShiftS s1 {p1, *ap.begin(), Point::meas(*p1, **ap.begin())};
		ShiftS s2 {p2, *ap.rbegin(), Point::meas(*p2, **ap.rbegin())};
		allshifts.push_back(s1);
		allshifts.push_back(s2);
	}
	//leave only point on outer grid contour which are non-significant for cont
	for (auto& s: allshifts){
		for (auto p: proc.cp){
			if (Point::meas(*s.from, *p) < geps*geps){
				s.from = NULL;
				break;
			}
		}
	}
	allshifts.remove_if([](ShiftS s){ return s.from==NULL; });
	//sort allshifts
	allshifts.sort();
	//make shifts
	auto point_cell = Connectivity::VertexCell(grid.vcells);
	aa::enumerate_ids_pvec(grid.vvert);
	for (auto& s: allshifts) if (s.from != NULL){
		//try to shift
		Point bu(*s.from);
		s.from->set(*s.to);
		//check if all cells are still ok
		auto& ps = point_cell[s.from->id];
		bool ok = true;
		for (auto i=0; i<ps.size(); ++i){
			auto& cc = grid.vcells[ps.cind[i]]->edges;
			if (std::get<0>(Contour::Finder::SelfCross(cc))){
				ok = false;
				break;
			}
		}
		if (ok){
			for (auto& s2: allshifts) if (&s2 != &s)
				if (s.to == s2.to || s.from == s2.from) s2.from = NULL;
		} else {
			s.from->set(bu);
		}
	}
}

namespace{
bool check_cells(const GridData& g, double& sumarea){
	//check cell-by-cell
	sumarea = 0;
	for (auto c: g.vcells){
		if (!Contour::IsContour(c->edges) || Contour::IsOpen(c->edges)) return false;
		if (std::get<0>(Contour::Finder::SelfCross(c->edges))) return false;
		double a = Contour::Area(c->edges);
		if (a <= 0) return false;
		sumarea += a;
	}
	return true;
}

bool check_edges_intersections(const GridData& g){
	auto bb = HM2D::BBox(g.vvert);
	BoundingBoxFinder bfinder(bb, bb.maxlen()/40);
	vector<BoundingBox> ebb(g.vedges.size());
	for (int i=0; i<g.vedges.size(); ++i){
		ebb[i] = BoundingBox(*g.vedges[i]->pfirst(), *g.vedges[i]->plast());
		bfinder.addentry(ebb[i]);
	}
	double ksieta[2];
	for (int i=0; i<g.vedges.size(); ++i)
	for (int isus: bfinder.suspects(ebb[i])) if (isus > i){
		auto& e0 = g.vedges[i];
		auto& e1 = g.vedges[isus];
		if (e0->pfirst() == e1->pfirst() || e0->pfirst() == e1->plast()) continue;
		if (e0->plast() == e1->pfirst() || e0->plast() == e1->plast()) continue;
		if (SectCross(*e0->pfirst(), *e0->plast(), *e1->pfirst(), *e1->plast(), ksieta)){
			return false;
		}
	}
	return true;
}
};

bool Algos::Check(const GridData& g){
	double sumarea = 0;
	if (!check_cells(g, sumarea)) return false;

	//check areas
	auto conts = HM2D::Contour::Assembler::GridBoundary(g);
	for (auto& c: conts) if (Contour::IsOpen(c)) return false;
	double tarea = 0;
	vector<EdgeData*> inners, outers;
	for (auto& c: conts){
		double ar = Contour::Area(c);
		int tt = 0;
		if (!(ar>0)) tt+=1;
		if (!Contour::CorrectlyDirectedEdge(c, 0)) tt+=2;
		if (!c[0]->has_left_cell()) tt += 4;
		switch (tt){
			case 0: case 5: case 6: case 3: ar = fabs(ar); break;
			case 4: case 1: case 2: case 7: ar = -fabs(ar); break;
		}
		tarea += ar;
		if (ar > 0) outers.push_back(&c);
		else inners.push_back(&c);
	}
	if (fabs(sumarea - tarea) > geps) return false;
	if (outers.size() == 0) return false;

	Contour::Tree sumtree;
	sumtree.add_contour(std::move(*outers[0]));
	for (int i=1; i<outers.size(); ++i){
		sumtree = Contour::Clip::Union(sumtree, *outers[i]);
	}
	for (int i=0; i<inners.size(); ++i){
		sumtree = Contour::Clip::Difference(sumtree, *inners[i]);
	}
	double tarea2 = sumtree.area();
	if (fabs(sumarea - tarea2) > geps) return false;

	if (fabs(sumarea - tarea) < geps*geps &&
	    fabs(sumarea - tarea2) < geps*geps) return true;
	
	//check edge-by-edge because inconcistency between areas
	//may occur as a result of numerical errors
	return check_edges_intersections(g);
}

namespace{

void remove_unused_prims(GridData& from){
	aa::constant_ids_pvec(from.vedges, 0);
	aa::constant_ids_pvec(from.vvert, 0);
	for (auto c: from.vcells)
	for (auto e: c->edges){
		e->id = 1;
		for (auto v: e->vertices) v->id = 1;
	}

	//remove unused edges
	//remove unused points
	aa::keep_by_id(from.vedges, 1);
	aa::keep_by_id(from.vvert, 1);
}

}

//merges coincident vertices, checks rotation
void Algos::Heal(GridData& from){
	//fix rotation
	for (int i=0; i<from.vcells.size(); ++i){
		auto& c = from.vcells[i];
		auto& ec = c->edges;
		if (Contour::Area(ec) < 0){ Contour::Reverse(ec); }
		//fix left-right cell connections
		for (int j=0; j<ec.size(); ++j){
			bool iscor = Contour::CorrectlyDirectedEdge(ec, j);
			if (iscor) ec[j]->left = c;
			else ec[j]->right = c;
		}
	}

	auto _less = [](const shared_ptr<Vertex>& p1,
			const shared_ptr<Vertex>& p2)->bool{
		return *p1 < *p2;
	};

	//--- merge congruent vertices
	VertexData av = from.vvert;
	std::sort(av.begin(), av.end(), _less);
	auto ur = std::unique(av.begin(), av.end(),
		[](const shared_ptr<Vertex>& a, const shared_ptr<Vertex>& b)
			{return *a == *b; });
	if (ur == av.end()) return;

	av.resize(ur - av.begin());
	aa::constant_ids_pvec(from.vvert, 0);
	aa::constant_ids_pvec(av, 1);

	for (auto& e: from.vedges)
	for (auto& v: e->vertices){
		if (v->id == 0) {
			v = *std::lower_bound(av.begin(), av.end(), v, _less);
			v->id = 2;
		}
	}

	//--- merge congruent edges
	//zero size edges and suspicios edges
	EdgeData zedges, susedges;
	for (auto& e: from.vedges)
	if (e->vertices[0]->id == 2 && e->vertices[1]->id == 2){
		if (e->vertices[0] == e->vertices[1]) zedges.push_back(e);
		else susedges.push_back(e);
	}

	//equal edges list
	typedef std::pair<Point*, Point*> Tppp;
	auto less_tppp = [](const Tppp* a, const Tppp* b)->bool{
		if (a->first == b->first) return a->second < b->second;
		else return a->first < b->first;
	};
	vector<Tppp> sus_sorted;
	for (auto e: susedges){
		auto pp = std::make_pair(e->first().get(), e->last().get());
		if (pp.first > pp.second) std::swap(pp.first, pp.second);
		sus_sorted.push_back(pp);
	}
	vector<Tppp*> psus_sorted;
	for (auto& it: sus_sorted) psus_sorted.push_back(&it);
	std::sort(psus_sorted.begin(), psus_sorted.end(), less_tppp);

	EdgeData change1, change2;
	auto it = psus_sorted.begin();
	while (it!= psus_sorted.end() && it != std::prev(psus_sorted.end())){
		auto v1 = *it++;
		auto v2 = *it;
		if (less_tppp(v1, v2)) continue;
		int ind1 = v1 - &sus_sorted[0];
		int ind2 = v2 - &sus_sorted[0];
		change1.push_back(susedges[ind1]);
		change2.push_back(susedges[ind2]);
		++it;
	}

	//edge->cell tables
	for (int i=0; i<change1.size(); ++i){
		auto& efrom = change2[i];
		auto& eto = change1[i];
		bool has_same_dir = efrom->first() == eto->first();
		if (has_same_dir){
			if (eto->no_left_cell()) eto->left = efrom->left;
			if (eto->no_right_cell()) eto->right = efrom->right;
		} else {
			if (eto->no_left_cell()) eto->left = efrom->right;
			if (eto->no_right_cell()) eto->right = efrom->left;
		}
	}

	//set edges ids
	aa::constant_ids_pvec(from.vedges, -1);
	aa::constant_ids_pvec(zedges, -2);
	aa::enumerate_ids_pvec(change2);
	
	//cell->edges table
	for (auto c: from.vcells){
		auto it = c->edges.begin();
		while (it != c->edges.end()){
			if ((*it)->id >= 0) (*it) = change1[(*it)->id];
			if ((*it)->id == -2) it = c->edges.erase(it);
			else ++it;
		}
	}

	//removing zero sized cells
	auto cit = std::remove_if(from.vcells.begin(), from.vcells.end(),
		[](shared_ptr<Cell>& c){ return c->edges.size() < 2; });
	from.vcells.resize(cit-from.vcells.begin());

	remove_unused_prims(from);
}

void Algos::RestoreFromCells(GridData& g){
	auto r1 = std::remove_if(g.vcells.begin(), g.vcells.end(),
			[](shared_ptr<Cell> c){ return c==nullptr; });
	g.vcells.resize(r1 - g.vcells.begin());

	//add primitives which are not in grid
	auto ae = AllEdges(g.vcells);
	auto av = AllVertices(ae);

	aa::constant_ids_pvec(ae, 0);
	aa::constant_ids_pvec(g.vedges, 1);
	aa::constant_ids_pvec(av, 0);
	aa::constant_ids_pvec(g.vvert, 1);

	std::copy_if(ae.begin(), ae.end(), std::back_inserter(g.vedges),
		[](shared_ptr<Edge> e){ return e->id == 0; });
	std::copy_if(av.begin(), av.end(), std::back_inserter(g.vvert),
		[](shared_ptr<Vertex> v){ return v->id == 0; });

	//remove unused primitives
	aa::constant_ids_pvec(ae, 0);
	aa::constant_ids_pvec(av, 0);
	aa::keep_by_id(g.vedges, 0);
	aa::keep_by_id(g.vvert, 0);
}

void Algos::RemoveShortEdges(GridData& grid, double ref_len){
	//1) calculate charachteristic cell sizes
	vector<double> csz(grid.vcells.size());
	for (int i=0; i<grid.vcells.size(); ++i){
		csz[i] = ref_len*sqrt(Contour::Area(grid.vcells[i]->edges));
	}
	//2) boundary points indicies
	aa::constant_ids_pvec(grid.vvert, 0);
	for (auto e: grid.vedges) if (e->is_boundary()){
		e->vertices[0]->id = 1;
		e->vertices[1]->id = 1;
	}
	//3) collapse short edges
	aa::enumerate_ids_pvec(grid.vcells);
	vector<int> badedges;
	for (int i=0; i<grid.vedges.size(); ++i) if (grid.vedges[i]->is_inner()){
		auto& e = grid.vedges[i];
		double elen = e->length();
		if (elen > csz[e->left.lock()->id]) continue;
		if (elen > csz[e->right.lock()->id]) continue;

		auto &p1 = e->vertices[0], &p2 = e->vertices[1];
		bool isbnd1 = p1->id == 1;
		bool isbnd2 = p2->id == 1;
		if (!isbnd1 && !isbnd2){
			double x = (p1->x + p2->x)/2.0, y = (p1->y + p2->y)/2.0;
			p1->set(x, y); p2->set(x, y);
		}
		else if (isbnd1 && !isbnd2) p2->set(*p1);
		else if (!isbnd1 && isbnd2) p1->set(*p2);
		badedges.push_back(i);
	}
	
	//4) merge points
	Heal(grid);
}

void Algos::NoConcaveCells(GridData& grid, double angle0, bool ignore_bnd){
	angle0 = angle0/180*M_PI;
	for (int i=0; i<grid.vcells.size(); ++i){
		if (grid.vcells[i]->edges.size()<4) continue;
		auto op = Contour::OrderedPoints(grid.vcells[i]->edges);
		auto opprev = op.end()[-2];
		for (int j=0; j<op.size()-1; ++j){
			if (ignore_bnd && grid.vcells[i]->edges[j]->is_boundary())
				continue;
			double ang = Angle(*opprev, *op[j], *op[j+1]);
			if (ISEQGREATER(ang, angle0)){
				if (SplitCell(grid, i, j)){
					--i;
					goto NEXT;
				}
			}
			opprev = op[j];
		}
	NEXT:
		continue;
	}
}
