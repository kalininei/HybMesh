#include "modgrid.hpp"
#include "modcont.hpp"
#include "healgrid.hpp"
#include <stack>
#include "assemble2d.hpp"
#include "contabs2d.hpp"
#include "finder2d.hpp"
#include "treverter2d.hpp"
#include "infogrid.hpp"
#include "debug_grid2d.hpp"
using namespace HM2D;
using namespace HM2D::Grid;

void Algos::UniqueRearrange(GridData& from){
	auto _vless = [](const shared_ptr<Vertex>& a, const shared_ptr<Vertex>& b){
				return *a < *b; };
	auto _eless = [](const shared_ptr<Edge>& a, const shared_ptr<Edge>& b){
				return a->pfirst()->id < b->pfirst()->id; };
	auto _cless = [](const shared_ptr<Cell>& a, const shared_ptr<Cell>& b){
				return a->edges[0]->id < b->edges[0]->id; };
	//vertices
	std::sort(from.vvert.begin(), from.vvert.end(), _vless);
	//edges
	aa::enumerate_ids_pvec(from.vvert);
	for (auto& e: from.vedges){
		if (e->first()->id > e->last()->id) e->reverse();
	}
	std::sort(from.vedges.begin(), from.vedges.end(), _eless);
	//cells
	aa::enumerate_ids_pvec(from.vedges);
	for (auto& c: from.vcells){
		int imin=0;
		int vmin=c->edges[0]->id;
		for (int i=1; i<c->edges.size(); ++i){
			if (c->edges[i]->id < vmin){
				imin = i; vmin = c->edges[i]->id;
			}
		}
		std::rotate(c->edges.begin(), c->edges.begin()+imin, c->edges.end());
	}
	std::sort(from.vcells.begin(), from.vcells.end(), _cless);
}

//adds non-repeating (by pointer) vertices, edges and cells
void Algos::ShallowAdd(const GridData& from, GridData& to){
	aa::constant_ids_pvec(from.vcells, 1);
	aa::constant_ids_pvec(from.vedges, 1);
	aa::constant_ids_pvec(from.vvert, 1);
	aa::constant_ids_pvec(to.vcells, 0);
	aa::constant_ids_pvec(to.vedges, 0);
	aa::constant_ids_pvec(to.vvert, 0);

	std::copy_if(from.vcells.begin(), from.vcells.end(),
			std::back_inserter(to.vcells),
			[](const shared_ptr<Cell>& c){
				return c->id != 0;
			});
	std::copy_if(from.vedges.begin(), from.vedges.end(),
			std::back_inserter(to.vedges),
			[](const shared_ptr<Edge>& c){
				return c->id != 0;
			});
	std::copy_if(from.vvert.begin(), from.vvert.end(),
			std::back_inserter(to.vvert),
			[](const shared_ptr<Vertex>& c){
				return c->id != 0;
			});
}

void Algos::RemoveCells(GridData& grid, vector<int> icells){
	std::sort(icells.begin(), icells.end());
	icells.resize(std::unique(icells.begin(), icells.end()) - icells.begin());
	//remove not used cells
	CellData notused;
	for (int i: icells) notused.push_back(grid.vcells[i]);
	aa::constant_ids_pvec(grid.vcells, 0);
	aa::constant_ids_pvec(notused, 1);
	aa::remove_by_id(grid.vcells, 1);
	//edge->cell connections
	EdgeData enotused;
	for (auto c: notused)
	for (auto e: c->edges){
		if (e->left.lock() == c) e->left.reset();
		if (e->right.lock() == c) e->right.reset();
		if (e->no_left_cell() && e->no_right_cell())
			enotused.push_back(e);
	}
	//remove unused edges
	if (enotused.size() == 0) return;
	aa::constant_ids_pvec(grid.vedges, 0);
	aa::constant_ids_pvec(enotused, 1);
	aa::remove_by_id(grid.vedges, 1);

	//remove unused vertices
	aa::constant_ids_pvec(AllVertices(enotused), 1);
	for (auto c: grid.vcells)
	for (auto e: c->edges){
		e->vertices[0]->id = 0;
		e->vertices[1]->id = 0;
	}
	aa::remove_by_id(grid.vvert, 1);
}
void Algos::RemoveCellsById(GridData& grid, int id){
	std::vector<int> badcells;
	for (int i=0; i<grid.vcells.size(); ++i)
	if (grid.vcells[i]->id == id) badcells.push_back(i);
	RemoveCells(grid, badcells);
}

void Algos::RemoveCells(GridData& grid, const Contour::Tree& domain, int what){
	CellData goodcells;
	if (what==INSIDE) goodcells = ExtractCells(grid, domain, OUTSIDE);
	if (what==OUTSIDE) goodcells = ExtractCells(grid, domain, INSIDE);
	std::swap(grid.vcells, goodcells);
	RestoreFromCells(grid);
}

void Algos::MergeBoundaries(const GridData& from, GridData& to){
	//== force 'from' boundary nodes lying on 'to' boundary to 'to'
	//assembling boundaries
	auto bvfrom = AllVertices(ECol::Assembler::GridBoundary(from));
	auto bedto = ECol::Assembler::GridBoundary(to);
	auto bbox = HM2D::BBox(bedto);
	//quick search for 'to' boundary edges
	BoundingBoxFinder tofinder(bbox, bbox.maxlen()/30.);
	for (auto e: bedto){
		tofinder.addentry(BoundingBox(*e->pfirst(), *e->plast()));
	}
	for (auto v: bvfrom){
		double ksi;
		for (auto ecand: tofinder.suspects(*v)){
			auto e = bedto[ecand].get();
			//if 'from' vertex equals 'to' vertex ignore it.
			//It will be splitted in MergeTo procedure
			if (*v == *e->pfirst() || *v == *e->plast()) break;
			isOnSection(*v, *e->pfirst(), *e->plast(), ksi);
			//if 'from' vertex lie in the middle of 'to' edge
			//split 'to' edge by this vertex
			if (ISIN_NN(ksi, 0, 1)){
				aa::enumerate_ids_pvec(to.vedges);
				SplitEdge(to, e->id, {*v}, true);
				bedto.push_back(to.vedges.back());
				tofinder.addentry(BoundingBox(*to.vedges.back()->pfirst(),
				                              *to.vedges.back()->plast()));
			}
		}
	}
	MergeTo(from, to);
}

void Algos::MergeTo(const GridData& from, GridData& to){
	auto bedfrom = ECol::Assembler::GridBoundary(from);
	auto bedto = ECol::Assembler::GridBoundary(to);
	auto bvfrom = AllVertices(bedfrom);
	auto bvto = AllVertices(bedto);

	//equal vertices
	auto to_vertedge = Connectivity::VertexEdge(to.vedges, bvto);
	Finder::VertexMatch finder(bvfrom);
	aa::constant_ids_pvec(bvfrom, 0);
	aa::constant_ids_pvec(bvto, 0);
	//change all equal vertices in 'to' grid to their 'from' pointers.
	for (auto& ve: to_vertedge){
		auto fndres = finder.find(*ve.v);
		if (!fndres) continue;
		for (auto ei: ve.eind){
			auto toe = to.vedges[ei];
			if (toe->vertices[0] == ve.v) toe->vertices[0] = fndres;
			if (toe->vertices[1] == ve.v) toe->vertices[1] = fndres;
		}
		fndres->id = 1;
	}

	//!!! 'from' should not contain duplicated boundary nodes
	assert([&](){
		for (auto& v: bvfrom){
			if (v != finder.find(*v)) return false;
		}
		return true;
	}());

	//suspicious edges
	EdgeData to_susp, from_susp; 
	for (auto e: bedto){
		if (e->first()->id == 1 && e->last()->id == 1) to_susp.push_back(e);
	}
	for (auto e: bedfrom){
		if (e->first()->id == 1 && e->last()->id == 1) from_susp.push_back(e);
	}

	//find out zero sized edges
	EdgeData zedges;
	auto check_zedge = [&zedges](shared_ptr<Edge> e)->bool{
			if (e->vertices[0] == e->vertices[1]){
				zedges.push_back(e);
				return true;
			}
			return false;
		};
	auto z1 = std::remove_if(to_susp.begin(), to_susp.end(), check_zedge);
	auto z2 = std::remove_if(from_susp.begin(), from_susp.end(), check_zedge);
	to_susp.resize(z1 - to_susp.begin());
	from_susp.resize(z2 - from_susp.begin());

	//temporary revert to_susp and from_susp for easy search
	//using copied vectors because order of edges in to_* will be changed
	//and reverter destructor would not work correctly otherwise.
	auto to_susp2 = to_susp, from_susp2 = from_susp;
	Contour::R::LeftCells rev1(to_susp2);
	Contour::R::LeftCells rev2(from_susp2);

	//sort susp arrays. After sorting they should match one by one
	auto esortkey0 = [](shared_ptr<Edge> a, shared_ptr<Edge> b){
		return a->vertices[0] < b->vertices[0];
	};
	auto esortkey1 = [](shared_ptr<Edge> a, shared_ptr<Edge> b){
		return a->vertices[1] < b->vertices[1];
	};
	std::sort(to_susp.begin(), to_susp.end(), esortkey0);
	std::sort(from_susp.begin(), from_susp.end(), esortkey1);
	//force matching
	for (int i=0; i<to_susp.size(); ++i){
		assert(i<from_susp.size());
		auto& eto = to_susp[i];
		int j = i;
		while (j<from_susp.size()){
			auto& efrom = from_susp[j];
			if (efrom->pfirst() == eto->plast() &&
			    efrom->plast() == eto->pfirst()){
				break;
			}
			++j;
		}
		if (j==from_susp.size()){
			to_susp.erase(to_susp.begin() + i);
			--i;
		} else if (j!=i){
			from_susp.erase(from_susp.begin()+i, from_susp.begin()+j);
		}
	}
	if (from_susp.size()>to_susp.size()) from_susp.resize(to_susp.size());
	assert(to_susp.size() == from_susp.size());

	//remove zero sized edges
	for (auto e: zedges){
		auto c = (e->has_left_cell()) ? e->left.lock()
		                              : e->right.lock();
		for (auto it = c->edges.begin(); it != c->edges.end(); ++it)
			if (*it == e) { c->edges.erase(it); break; }
	}
	//split edges
	for (int i=0; i<from_susp.size(); ++i){
		auto c1 = to_susp[i]->left.lock();
		auto c2 = from_susp[i]->left.lock();
		//cell->edge
		for (int j=0; j<c1->edges.size(); ++j){
			if (c1->edges[j] == to_susp[i]){
				c1->edges[j] = from_susp[i];
				break;
			}
		}
		//edge->cell
		from_susp[i]->right = c1;
	}

	//add and reassemble
	to.vcells.insert(to.vcells.end(), from.vcells.begin(), from.vcells.end());
	RestoreFromCells(to);
}

namespace{

vector<bool> mark_segments(const EdgeData& seg, double angle){
	//seg contour was really directed
	vector<bool> ret(seg.size());
	VertexData keep_pts = AllVertices(ECol::Algos::Simplified(seg, angle));
	//edges starting with keep_pts -> id = 1; else id = 0;
	for (int i=0; i<seg.size(); ++i){
		auto e = seg[i];
		if (Finder::Contains(keep_pts, e->pfirst())) ret[i]=true;
		else ret[i]=false;
	}
	return ret;
}

//marks all unused boundaries with -1;
void simplify_boundary(Cell* cell, double angle){
	HM2D::Contour::R::ReallyDirect rd(cell->edges);
	//1) connect boundary edges
	vector<HM2D::EdgeData> be;
	for (auto e: cell->edges) if (e->is_boundary()){
		if (be.size() == 0 ||
		    be.back().back()->plast() != e->pfirst()){
			be.emplace_back();
		}
		be.back().push_back(e);
	}
	//check first last connection
	if (be.size()>1 && be[0][0]->pfirst() == be.back().back()->plast()){
		std::copy(be[0].begin(), be[0].end(), std::back_inserter(
				be.back()));
		std::swap(be[0], be.back());
		be.resize(be.size()-1);
	}
	//2) do simplifications
	vector<vector<bool>> keepbe;
	for (auto& it: be) keepbe.push_back(mark_segments(it, angle));

	//3) mark unused with id=-1
	aa::enumerate_ids_pvec(cell->edges);
	for (int i=0; i<be.size(); ++i)
	for (int j=0; j<be[i].size(); ++j){
		if (!keepbe[i][j]) cell->edges[be[i][j]->id]->id = -1;
	}
	//4) check resulting contour
	VertexData tp;
	for (auto& e: cell->edges) if (e->id != -1){
		tp.push_back(e->first());
	}
	EdgeData newc = Contour::Assembler::Contour1(tp, true);
	if (newc.size() < 3 ||
	    std::get<0>(Contour::Finder::SelfCross(newc)) ||
	    Contour::Area(newc) < 0){
		//simplification is invalid
		aa::constant_ids_pvec(cell->edges, 1);
		return;
	} else {
		//move vertices of valid edges
		tp.push_back(tp[0]);
		int iv=0;
		for (auto& e: cell->edges) if (e->id !=-1){
			e->vertices[1] = tp[++iv];
		}
	}
}

};

void Algos::SimplifyBoundary(GridData& grid, double angle){
	//list of cells which have more then two boundary edges
	vector<Cell*> bcells;
	vector<int> used(grid.vcells.size(), 0);
	aa::enumerate_ids_pvec(grid.vcells);
	for (auto e: grid.vedges) if (e->is_boundary()){
		auto c = (e->has_left_cell()) ? e->left.lock()
		                              : e->right.lock();
		++used[c->id];
		if (used[c->id] == 2) bcells.push_back(c.get());
	}
	//do simplifications
	for (auto c: bcells){
		simplify_boundary(c, angle);
		aa::remove_by_id(c->edges, -1);
	}
	//remove unused primitives
	RestoreFromCells(grid);
}


namespace{
//here we use cont which was already really directed
//in a counterclockwise direction
int find_worst_vertex(const EdgeData& cont){
	vector<double> angles(cont.size());
	for (int i=0; i<cont.size(); ++i){
		int ip = (i==0) ? cont.size()-1 : i-1;
		angles[i] = Angle(*cont[ip]->first(), *cont[i]->first(), *cont[i]->last());
	}
	return std::max_element(angles.begin(), angles.end()) - angles.begin();
}
bool check_connection(const EdgeData& cont, int n1, int n2){
	//if (cont.size()<5) return true;
	Point* p1 = cont[n1]->first().get();
	Point* p2 = cont[n2]->first().get();
	double ksi[2];
	for (int i=0; i<cont.size(); ++i){
		int ip = (i==cont.size()-1) ? 0 : i+1;
		if (i==n1 || i==n2) continue;
		if (ip==n1 || ip==n2) continue;
		if (SectCross(*p1, *p2, *cont[i]->vertices[0],
		              *cont[i]->vertices[1], ksi)){
			return false;
		}
	}
	return true;
}
int find_best_connection(const EdgeData& cont, int n){
	double ret=-1;
	double nrm = 1e54;
	int np = (n == 0) ? cont.size()-1 : n-1;
	int nn = (n == cont.size()-1) ? 0 : n+1;
	//a0 is the best angle.
	double a0 = Angle(*cont[np]->first(), *cont[n]->first(), *cont[nn]->first())/2.0;
	for (int i=n+2; i<=cont.size()+n-2; ++i){
		int ii = i % cont.size();
		if (!check_connection(cont, n, ii)) continue;
		//check if connection lies within contour
		double ian = Angle(*cont[np]->first(), *cont[n]->first(), *cont[ii]->first());
		if (ISEQ(ian, 0) || ISEQGREATER(ian, 2*a0)) continue;
		//calculate difference from a0. The lower the better.
		double diff = fabs(a0-ian);
		if (diff<nrm){nrm=diff; ret=ii;}
	}
	return ret;
}
}
bool Algos::SplitCell(GridData& grid, int icell, int lnode1, int lnode2){
	shared_ptr<Cell> oldcell(grid.vcells[icell]);
	EdgeData gc = oldcell->edges;
	Contour::R::ReallyDirect rev(gc);

	//calculate lnodes
	if (lnode1 < 0) lnode1 = find_worst_vertex(gc);
	if (lnode2 < 0) lnode2 = find_best_connection(gc, lnode1);
	else if (!check_connection(gc, lnode1, lnode2)){
		lnode2 = -1;
	}
	if (lnode1 > lnode2) std::swap(lnode1, lnode2);
	//return false on bad lnode options
	if (lnode1 < 0 || lnode1 == lnode2 || lnode1>=gc.size() ||
	    lnode2>=gc.size() || lnode2-lnode1==1 || lnode2-lnode1==gc.size()-1) return false;
	
	//construct objects
	shared_ptr<Cell> newcell(new Cell);
	shared_ptr<Edge> newedge(new Edge(gc[lnode1]->first(), gc[lnode2]->first()));
	newedge->left = oldcell;
	newedge->right = newcell;
	oldcell->edges.erase(oldcell->edges.begin()+lnode1,
	                     oldcell->edges.begin()+lnode2);
	oldcell->edges.insert(oldcell->edges.begin()+lnode1, newedge);
	newcell->edges = vector<shared_ptr<Edge>>(gc.begin()+lnode1, gc.begin()+lnode2);
	for (auto& e: newcell->edges) e->left = newcell; 
	newcell->edges.push_back(newedge);
	
	//add to return
	grid.vcells.push_back(newcell);
	grid.vedges.push_back(newedge);
	return true;
}

bool Algos::SplitEdge(GridData& grid, int iedge, const vector<Point>& apoints, bool force){
	assert(iedge < grid.vedges.size());
	auto ed = grid.vedges[iedge];
	auto lc = ed->left.lock();
	auto rc = ed->right.lock();
	int lcloc = (!lc) ? -1 : std::find(lc->edges.begin(), lc->edges.end(), ed) -
	                         lc->edges.begin();
	int rcloc = (!rc) ? -1 : std::find(rc->edges.begin(), rc->edges.end(), ed) -
	                         rc->edges.begin();

	//build new vertices
	VertexData vert(apoints.size() + 2);
	vert[0] = ed->vertices[0];
	vert.back() = ed->vertices[1];
	for (int i=0; i<apoints.size(); ++i){
		vert[i+1].reset(new Vertex(apoints[i]));
	}
	//build new edges
	EdgeData newedges(apoints.size()+1);
	for (int i=0; i<apoints.size()+1; ++i){
		newedges[i].reset(new Edge(*ed));
		newedges[i]->vertices[0] = vert[i];
		newedges[i]->vertices[1] = vert[i+1];
	}
	//add edges to cells connectivity
	EdgeData nl, nr;
	if (lc){
		nl = lc->edges;
		nl[lcloc] = newedges[0];
		nl.insert(nl.begin()+lcloc+1, newedges.begin()+1, newedges.end());
	}
	if (rc){
		nr = rc->edges;
		nr[rcloc] = newedges[0];
		nr.insert(nr.begin()+rcloc, newedges.rbegin(), newedges.rend()-1);
	}

	//check for no cross
	if (!force){
		if (lc && std::get<0>(Contour::Finder::SelfCross(nl))) return false;
		if (rc && std::get<0>(Contour::Finder::SelfCross(nr))) return false;
	}

	//write calculated data to grid
	grid.vvert.insert(grid.vvert.end(), vert.begin()+1, vert.end()-1);
	grid.vedges.insert(grid.vedges.end(), newedges.begin()+1, newedges.end());
	ed->vertices[1] = vert[1];
	if (lc){
		nl[lcloc] = ed;
		lc->edges = nl;
	}
	if (rc){
		nr[rcloc+newedges.size()-1] = ed;
		rc->edges = nr;
	}

	return true;
}

//all cells with dimension higher than maxdim will be splitted
void Algos::CutCellDims(GridData& grid, int maxdim){
	assert(maxdim >= 3);

	for (int i=0; i<grid.vcells.size(); ++i){
		if (grid.vcells[i]->edges.size()>maxdim){
			if (SplitCell(grid, i)){
				//analyze i-th cell once again
				--i;
			}
		}
	}

}

void Algos::RemoveEdges(GridData& grid, vector<int> iedges){
	std::sort(iedges.begin(), iedges.end());
	iedges.resize(std::unique(iedges.begin(), iedges.end()) - iedges.begin());
	std::vector<int> badcells;
	aa::enumerate_ids_pvec(grid.vcells);
	for (auto ie: iedges){
		auto& e = grid.vedges[ie];
		auto cl = e->left.lock();
		auto cr = e->right.lock();
		if (cl==nullptr) badcells.push_back(cr->id);
		else if (cr==nullptr) badcells.push_back(cl->id);
		else{
			badcells.push_back(cr->id);
			int iloc_left = aa::shpvec_ifind(cl->edges, e.get());
			int iloc_right = aa::shpvec_ifind(cr->edges, e.get());
			std::rotate(cl->edges.begin(), cl->edges.begin()+iloc_left+1, cl->edges.end());
			std::rotate(cr->edges.begin(), cr->edges.begin()+iloc_right+1, cr->edges.end());
			cl->edges.resize(cl->edges.size()-1);
			cl->edges.insert(cl->edges.end(), cr->edges.begin(), cr->edges.end()-1);
			e->left.reset();
			for (int i=0; i<cl->edges.size(); ++i){
				if (cl->edges[i]->left.lock()==cr) cl->edges[i]->left = cl;
				else if (cl->edges[i]->right.lock()==cr) cl->edges[i]->right = cl;
			}
		}
	}
	RemoveCells(grid, badcells);
}
