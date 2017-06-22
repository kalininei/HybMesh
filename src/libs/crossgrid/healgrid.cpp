#include "healgrid.hpp"
#include "modcont.hpp"
#include "modgrid.hpp"
#include "treverter2d.hpp"
#include "assemble2d.hpp"
#include "contabs2d.hpp"
#include "finder2d.hpp"
#include "clipdomain.hpp"
#include "nodes_compare.h"

using namespace HM2D;
using namespace HM2D::Grid;

namespace{
bool check_cells(const GridData& g, double& sumarea){
	//check cell-by-cell
	sumarea = 0;
	for (auto c: g.vcells){
		if (!Contour::IsContour(c->edges) || Contour::IsOpen(c->edges)){
			return false;
		}
		if (std::get<0>(Contour::Finder::SelfCross(c->edges))){
			return false;
		}
		double a = Contour::Area(c->edges);
		if (a <= 0){
			return false;
		}
		sumarea += a;
	}
	return true;
}

bool check_edges_intersections(const GridData& g){
	auto bb = HM2D::BBox(g.vvert);
	BoundingBoxFinder bfinder(bb, bb.maxlen()/50);
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
	if (outers.size() < 2 && fabs(sumarea - tarea) < geps*geps) return true;

	//we can not relay only on sum-cell-areas == sum-bound-areas check
	//because grid can contain independant overlapping parts.
	//To treat the latter case we have to explicitly calculate edge intersections.
	//(... maybe check for bound intersections will be enough for most cases?)
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
		if (Contour::Area(ec) < 0){
			Contour::Algos::Reverse(ec);
		}
		//fix left-right cell connections
		for (int j=0; j<ec.size(); ++j){
			bool iscor = Contour::CorrectlyDirectedEdge(ec, j);
			if (iscor){
				if (ec[j]->right.lock() == c) std::swap(ec[j]->right, ec[j]->left);
				else ec[j]->left = c;
			} else {
				if (ec[j]->left.lock() == c) std::swap(ec[j]->right, ec[j]->left);
				else ec[j]->right = c;
			}
		}
	}

	auto _less = [](const shared_ptr<Vertex>& p1,
			const shared_ptr<Vertex>& p2)->bool{
		return *p1 < *p2;
	};

	//--- merge congruent vertices
	auto nset = point2_set_pointers(from.vvert);
	vector<std::pair<int, int>> from_to = nset.verbose_unique();
	if (from_to.size() == 0) return;

	aa::constant_ids_pvec(from.vvert, -1);
	for (auto it: from_to){
		from.vvert[it.first]->id = it.second;
		from.vvert[it.second]->id = -2;
	}

	for (auto& e: from.vedges)
	for (auto& v: e->vertices) 
		if (v->id > -1) v = from.vvert[v->id];

	//--- merge congruent edges
	//zero size edges and suspicios edges
	EdgeData zedges, susedges;
	for (auto& e: from.vedges)
	if (e->vertices[0]->id == -2 && e->vertices[1]->id == -2){
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
		if (ind1 > ind2) std::swap(ind1, ind2);
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
