#include "modgrid.hpp"
#include "algos.hpp"
#include "healgrid.hpp"
#include <stack>
#include "cont_assembler.hpp"
#include "contabs2d.hpp"
#include "finder2d.hpp"
#include "treverter2d.hpp"
using namespace HM2D;
using namespace HM2D::Grid;

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

namespace{
//returns list of fully outer b indicies
std::vector<int> filter_by_bbox(const vector<BoundingBox>& b, const BoundingBox& bbox){
	std::vector<int> ret;
	for (int i=0; i<b.size(); ++i){
		if (!bbox.has_common_points(b[i])) ret.push_back(i);
	}
	return ret;
}
vector<int> fill_group(int first, vector<bool>& used, int nx, int ny){
	std::array<int, 4> adj;
	auto fill_adj = [&nx, &ny, &adj](int i){
		int ix = i % nx;
		int iy = i / nx;
		adj[0] = (ix != 0) ? (iy*nx+ix-1) : -1;
		adj[1] = (ix != nx-1) ? (iy*nx+ix+1) : -1;
		adj[2] = (iy != 0) ? ((iy-1)*nx+ix) : -1;
		adj[3] = (iy != ny-1) ? ((iy+1)*nx+ix) : -1;
	};

	vector<int> ret(1, first); used[first] = true;
	int uu=0;
	while (uu<ret.size()){
		fill_adj(ret[uu]);
		for (auto a: adj) if (a != -1 && !used[a]){
			ret.push_back(a);
			used[a] = true;
		}
		++uu;
	}

	return ret;
}

vector<vector<int>> squares_group(const BoundingBoxFinder& fn, const vector<int>& tp){
	vector<bool> used(tp.size(), false);
	vector<vector<int>> ret;
	for (int i=0; i<tp.size(); ++i) if (tp[i] != -1) used[i] = true;

	while (1){
		int first = std::find(used.begin(), used.end(), false)-used.begin();
		if (first >= used.size()) break;
		ret.push_back(fill_group(first, used, fn.nx(), fn.ny()));
	}

	return ret;
}

vector<int> define_squares_positions(const BoundingBoxFinder& bf,
		const Contour::Tree& tree, const EdgeData& tree_edges){
	//0 - inactive, 1 - inside, 2 - outside, 3 - undefined
	vector<int> ret(bf.nsqr(), -1);

	//set undefined to those which contains tree edges
	for (auto e: tree_edges){
		for (int i: bf.suspects(BoundingBox(*e->first(), *e->last()))){
			ret[i] = 3;
		}
	}

	//set inactive to those which do not contain any grid cells
	for (int i=0; i<bf.nsqr(); ++i) if (bf.sqr_entries(i).size() == 0){
		ret[i] = 0;
	}

	//group all left squares
	vector<vector<int>> groups = squares_group(bf, ret);

	//get feature for each group
	vector<Point> gpoints;
	for (auto& g: groups){
		gpoints.push_back(bf.sqr_center(g[0]));
	}
	vector<int> gf = Contour::Algos::SortOutPoints(tree, gpoints);
	for (int i=0; i<groups.size(); ++i){
		for (auto ib: groups[i]) ret[ib] = gf[i];
	}

	return ret;
}

}

CellData Algos::ExtractCells(const GridData& grid, const Contour::Tree& domain, int what){
	//all good cells will be stored in ret
	CellData ret;
	//tree simplification
	Contour::Tree nt = Contour::Algos::Simplified(domain);
	EdgeData nted = nt.alledges();
	//calculate bounding boxes for all cells
	vector<BoundingBox> gridbb(grid.vcells.size());
	for (int i=0; i<grid.vcells.size(); ++i){
		gridbb[i] = BBox(grid.vcells[i]->edges);
	}

	//leave only those cells which lie inside domain box
	auto bbox2 = HM2D::BBox(nted);
	std::vector<int> outercells = filter_by_bbox(gridbb, bbox2);
	aa::constant_ids_pvec(grid.vcells, 0);
	for (auto i: outercells) grid.vcells[i]->id=1;
	CellData icells;
	std::copy_if(grid.vcells.begin(), grid.vcells.end(), std::back_inserter(icells),
			[](const shared_ptr<Cell>& c){ return c->id == 0; });
	VertexData ivert=AllVertices(icells);

	//add excluded cells to result if we want to exclude inner domain area
	if (what == OUTSIDE)
	std::copy_if(grid.vcells.begin(), grid.vcells.end(), std::back_inserter(ret),
			[](const shared_ptr<Cell>& c){ return c->id == 1; });

	//build finder
	BoundingBox bbox({HM2D::BBox(ivert), bbox2});
	BoundingBoxFinder bfinder(bbox, bbox.maxlen()/50);
	aa::enumerate_ids_pvec(grid.vcells);
	for (auto c: icells){
		bfinder.addentry(gridbb[c->id]);
	}

	//get square positions
	//0 - inactive, 1 - inside, 2 - outside, 3 - undefined
	vector<int> sqrpos = define_squares_positions(bfinder, nt, nted);

	//fill cell ids
	for (int i=0; i<sqrpos.size(); ++i) if (sqrpos[i] == 1 || sqrpos[i] == 2){
		for (int k: bfinder.sqr_entries(i)) icells[k]->id = sqrpos[i];
	}
	for (int i=0; i<sqrpos.size(); ++i) if (sqrpos[i] == 3){
		for (int k: bfinder.sqr_entries(i)) icells[k]->id = 3;
	}

	//add good cells to result
	int goodid = (what == INSIDE) ? 1 : 2;
	std::copy_if(grid.vcells.begin(), grid.vcells.end(), std::back_inserter(ret),
			[&goodid](const shared_ptr<Cell>& c){ return c->id == goodid; });

	//leave only undefined cells
	aa::keep_by_id(icells, 3);
	ivert = AllVertices(icells);

	//Sorting vertices
	vector<Point> pivert(ivert.size());
	for (int i=0; i<ivert.size(); ++i) pivert[i].set(*ivert[i]);
	vector<int> srt = Contour::Algos::SortOutPoints(nt, pivert);
	for (int i=0; i<ivert.size(); ++i) ivert[i]->id = srt[i];

	//analyzing undefined cells
	CellData bndcells;
	int badid = (what == INSIDE) ? OUTSIDE : INSIDE;
	for (auto c: icells){
		bool allbnd = true;
		for (auto e: c->edges)
		for (auto v: e->vertices){
			if (v->id == badid) goto BADCELL;
			if (allbnd && v->id != BOUND) allbnd = false;
		}
		//if all points lie on boundary
		if (allbnd){ 
			bndcells.push_back(c);
			goto BADCELL;
		}
		ret.push_back(c);
	BADCELL:
		continue;
	}
	//analyzing cells with all points lying on boundary
	for (auto c: bndcells){
		Point p = Contour::InnerPoint(c->edges);
		int pos = nt.whereis(p);
		if (pos == what) ret.push_back(c);
	}

	return ret;
}

void Algos::RemoveCells(GridData& grid, const vector<int>& icells){
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

void Algos::RemoveCells(GridData& grid, const Contour::Tree& domain, int what){
	CellData goodcells;
	if (what==INSIDE) goodcells = ExtractCells(grid, domain, OUTSIDE);
	if (what==OUTSIDE) goodcells = ExtractCells(grid, domain, INSIDE);
	std::swap(grid.vcells, goodcells);
	RestoreFromCells(grid);
}

void Algos::MergeTo(const GridData& from, GridData& to){
	auto bedfrom = ECol::Assembler::GridBoundary(from);
	auto bedto = ECol::Assembler::GridBoundary(to);

	//equal vertices
	auto to_vertedge = Connectivity::VertexEdge(bedto);
	Finder::VertexMatch finder(AllVertices(bedfrom));
	aa::constant_ids_pvec(bedfrom, 0);
	aa::constant_ids_pvec(bedto, 0);
	for (auto& ve: to_vertedge){
		auto fndres = finder.find(*ve.v);
		if (!fndres) continue;
		for (auto ei: ve.eind){
			auto toe = bedto[ei];
			if (toe->vertices[0] == ve.v) toe->vertices[0] = fndres;
			if (toe->vertices[1] == ve.v) toe->vertices[1] = fndres;
		}
		fndres->id = 1;
	}

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
	assert([&](){
			if (to_susp.size() != from_susp.size()) return false;
			for (int i=0; i<to_susp.size(); ++i){
				if (to_susp[i]->first() != from_susp[i]->last()) return false;
				if (to_susp[i]->last() != from_susp[i]->first()) return false;
			}
			return true;
		}());

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

void Algos::SimplifyBoundary(GridData& grid, double angle){
	//list of cells which have more then two boundary edges
	vector<Cell*> bcells;
	vector<bool> added(grid.vcells.size(), false);
	aa::enumerate_ids_pvec(grid.vcells);
	for (auto e: grid.vedges) if (e->is_boundary()){
		auto c = (e->has_left_cell()) ? e->left.lock()
		                              : e->right.lock();
		if (added[c->id]) bcells.push_back(c.get());
		else added[c->id] = true;
	}

	//get sequence of edges to connect
	vector<EdgeData> toconnect;
	for (auto c: bcells){
		EdgeData bed;
		for (auto e: c->edges) if (e->is_boundary())
			bed.push_back(e);
		vector<EdgeData> sc = Contour::Assembler::SimpleContours(bed);
		for (auto bcont: sc) if (bcont.size() > 1){
			EdgeData simp = ECol::Algos::Simplified(bcont, angle);
			VertexData keep_pnt(AllVertices(simp));
			auto op = Contour::OrderedPoints(bcont);
			if (keep_pnt.size() == op.size()) continue;
			EdgeData* addto = 0;
			for (int i=0; i<op.size()-1; ++i){
				if (Finder::Contains(keep_pnt, op[i].get())){
					toconnect.emplace_back();
					addto = &toconnect.back();
				}
				addto->push_back(bcont[i]);
			}
		}
	}
	//connect edges
	aa::constant_ids_pvec(grid.vedges, 1);
	aa::constant_ids_pvec(grid.vvert, 1);
	for (auto& it: toconnect) if (it.size()>1){
		//remove all connect edges except first one
		Contour::R::LeftCells::Permanent(it);
		auto c = it[0]->left.lock();
		EdgeData olde = c->edges;
		EdgeData be(it.begin()+1, it.end());
		auto r = std::remove_if(c->edges.begin(), c->edges.end(),
			[&be](shared_ptr<Edge> e){ return Finder::Contains(be, e.get()) != nullptr; });
		c->edges.resize(r - c->edges.end());
		it[0]->vertices[1] = it.back()->vertices[1];
		//if cell has self crosses revert all changes back
		if (std::get<0>(Contour::Finder::SelfCross(c->edges))){
			it[0]->vertices[1] = it[1]->vertices[0];
			c->edges = olde;
		} else 
		//else mark all deleted points and edges with zeros
		for (int i=1; i<it.size(); ++i){
			it[i]->vertices[0]->id = 0;
			it[i]->id = 0;
		}
	}
	//remove unused
	aa::remove_by_id(grid.vvert, 0);
	aa::remove_by_id(grid.vedges, 0);
}
namespace{
//here we use cont which was already really directed
//in a counterclockwise direction
int find_worst_vertex(const EdgeData& cont){
	vector<double> angles(cont.size());
	for (int i=0; i<cont.size(); ++i){
		int ip = (i==0) ? cont.size()-1 : i-1;
		double a = Angle(*cont[ip]->first(), *cont[i]->first(), *cont[i]->last());
		angles[i] = fabs(M_PI-a);
	}
	return std::max_element(angles.begin(), angles.end()) - angles.begin();
}
bool check_connection(const EdgeData& cont, int n1, int n2){
	if (cont.size()<5) return true;
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
	double a0 = Angle(*cont[np]->first(), *cont[n]->first(), *cont[nn]->first())/2.0;
	for (int i=n+2; i<=cont.size()+n-2; ++i){
		int ii = i % cont.size();
		if (!check_connection(cont, n, ii)) continue;
		double ian = Angle(*cont[np]->first(), *cont[n]->first(), *cont[ii]->first());
		if (ian >= 2*a0) continue;
		double diff = fabs(a0-ian);
		if (diff<nrm){nrm=diff; ret=i;}
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
	newcell->edges.push_back(newedge);
	
	//add to return
	grid.vcells.push_back(newcell);
	grid.vedges.push_back(newedge);
	return true;
}

bool Algos::SplitEdge(GridData& grid, int iedge, const vector<Point>& apoints, bool force){
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
	for (int i=0; i<apoints.size(); ++i){
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

//all cells with dimensin higher than maxdim will be splitted
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
