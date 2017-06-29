#include "infogrid.hpp"
#include "contour_tree.hpp"
#include "modcont.hpp"
#include "finder2d.hpp"

using namespace HM2D;
namespace hg=HM2D::Grid;

double hg::Area(const GridData& grid){
	auto vt = Contour::Tree::GridBoundary01(grid);
	double ret = 0;
	for (auto it: vt) ret += it.area();
	return ret;
}

vector<double> hg::CellAreas(const GridData& grid){
	vector<double> ret(grid.vcells.size());
	for (int i=0; i<ret.size(); ++i){
		ret[i] = HM2D::Contour::Area(grid.vcells[i]->edges);
	}
	return ret;
}

//calculate skewness
vector<double> hg::Skewness(const GridData& grid){
	vector<double> ret(grid.vcells.size());
	for (int i=0; i<grid.vcells.size(); ++i){
		const Cell& c = *grid.vcells[i];
		int dim = c.edges.size();
		if (dim < 3){
			ret[i] = 1.0;  //very bad cell anyway
			continue;
		}
		vector<double> angles(dim);
		auto op = Contour::OrderedPoints(c.edges);
		for (int j=1; j<dim; ++j){
			const Point& p0 = *op[j-1];
			const Point& p1 = *op[j];
			const Point& p2 = *op[j+1];
			angles[j] = Angle(p0, p1, p2);
		}
		angles[0] = M_PI *(dim - 2) - 
			std::accumulate(angles.begin() + 1, angles.begin() + dim, 0.0);
		auto minmax = std::minmax_element(angles.begin(), angles.end());
		double minv = *minmax.first;
		double maxv = *minmax.second;
		double refan = M_PI * (dim - 2) / dim;
		ret[i] = std::max( (maxv-refan)/(M_PI-refan), (refan-minv)/refan );
		if (ret[i] > 1.0) ret[i] = 1.0;   //for non-convex cells
	}
	return ret;
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
		for (int i: bf.sqrs_by_bbox(BoundingBox(*e->first(), *e->last()))){
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
	vector<int> gf = Contour::Finder::SortOutPoints(tree, gpoints);
	for (auto& i: gf){
		if (i==INSIDE) i = 1;
		else if (i==OUTSIDE) i = 2;
		else i = 3;
	}
	for (int i=0; i<groups.size(); ++i){
		for (auto ib: groups[i]) ret[ib] = gf[i];
	}

	return ret;
}

}

HMCallback::FunctionWithCallback<Grid::TExtractCells> Grid::ExtractCells;

CellData Grid::TExtractCells::_run(const GridData& grid, const Contour::Tree& domain, int what){
	if (domain.roots().size() == 0){
		if (what == INSIDE || what == BOUND) return {};
		else return grid.vcells;
	}
	if (grid.vcells.size() == 0) return {};
	callback->step_after(20, "Box precessing", 8, 1);
	//all good cells will be stored in ret
	CellData ret;
	//tree simplification
	Contour::Tree nt = Contour::Algos::Simplified(domain);
	if (what != BOUND) nt.remove_detached();
	EdgeData nted = nt.alledges();
	//calculate bounding boxes for all cells
	vector<BoundingBox> gridbb(grid.vcells.size());
	for (int i=0; i<grid.vcells.size(); ++i){
		gridbb[i] = BBox(grid.vcells[i]->edges);
	}

	//leave only those cells which lie inside domain box
	callback->subprocess_step_after(1);
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

	//return if no cells to process
	if (icells.size() == 0) return ret;

	//build finder
	callback->subprocess_step_after(2);
	BoundingBox bbox({HM2D::BBox(ivert), bbox2});
	BoundingBoxFinder bfinder(bbox, bbox.maxlen()/30);
	aa::enumerate_ids_pvec(grid.vcells);
	for (auto c: icells){
		bfinder.addentry(gridbb[c->id]);
	}

	//get square positions
	//0 - inactive, 1 - inside, 2 - outside, 3 - undefined
	callback->subprocess_step_after(3);
	vector<int> sqrpos = define_squares_positions(bfinder, nt, nted);

	//fill cell ids
	callback->subprocess_step_after(1);
	for (int i=0; i<sqrpos.size(); ++i) if (sqrpos[i] == 1 || sqrpos[i] == 2){
		for (int k: bfinder.sqr_entries(i)) icells[k]->id = sqrpos[i];
	}
	for (int i=0; i<sqrpos.size(); ++i) if (sqrpos[i] == 3){
		for (int k: bfinder.sqr_entries(i)) icells[k]->id = 3;
	}

	if (what == BOUND){
		callback->step_after(80, "Matching cells");
		aa::keep_by_id(icells, 3);
		return extract_bound(icells, nted);
	}

	//add good cells to result
	int goodid = (what == INSIDE) ? 1 : 2;
	std::copy_if(icells.begin(), icells.end(), std::back_inserter(ret),
			[&goodid](const shared_ptr<Cell>& c){ return c->id == goodid; });

	//leave only undefined cells
	aa::keep_by_id(icells, 3);
	ivert = AllVertices(icells);

	//Sorting vertices
	callback->step_after(25, "Sorting points");
	vector<Point> pivert(ivert.size());
	for (int i=0; i<ivert.size(); ++i) pivert[i].set(*ivert[i]);
	vector<int> srt = Contour::Finder::SortOutPoints(nt, pivert);
	for (int i=0; i<ivert.size(); ++i) ivert[i]->id = srt[i];

	//analyzing undefined cells: if it contains bad point
	//it is undoubtedly bad, else add cell for the next check query (suspcells).
	callback->step_after(55, "Sorting cells", 10, 7);
	CellData suspcells;
	int badid = (what == INSIDE) ? OUTSIDE : INSIDE;
	for (auto c: icells){
		//vertices check
		for (auto e: c->edges)
		for (auto v: e->vertices){
			if (v->id == badid) goto CONTINUE1;
		}
		//inner point check
		{
			Point ip = Contour::InnerPoint(c->edges);
			if (domain.whereis(ip) == badid) goto CONTINUE1;
		}

		//add cell for future checks
		suspcells.push_back(c);
CONTINUE1:
		continue;
	}

	//even if all points are good, cell could be not fully good.
	//we do explicit intersection check here to be sure.
	callback->subprocess_step_after(3);
	if (suspcells.size() > 0){
		EdgeData suspedges = AllEdges(suspcells);
		aa::constant_ids_pvec(suspedges, 0);
		BoundingBox bb=HM2D::BBox(nted);
		BoundingBoxFinder domainfinder(bb, bb.maxlen()/30);
		for (auto& e: nted){
			domainfinder.addentry(BoundingBox(*e->pfirst(), *e->plast()));
		}
		double ksieta[2];
		for (auto& e: suspedges){
			//shrink edge to eliminate non-crossing intersection
			Point p1 = Point::Weigh(*e->pfirst(), *e->plast(), 2.*geps);
			Point p2 = Point::Weigh(*e->pfirst(), *e->plast(), 1.-2.*geps);
			for (auto isus: domainfinder.suspects(BoundingBox(p1, p2))){
				Point *n1 = nted[isus]->pfirst(), *n2 = nted[isus]->plast();
				//ignore parallel sections
				int w1 = LinePointWhereIs(p1, *n1, *n2);
				int w2 = LinePointWhereIs(p2, *n1, *n2);
				if ((w1 == 0 && w2 == 2) || (w2 == 0 && w1 == 2))
				if (SectCross(*n1, *n2, p1, p2, ksieta)){
					e->id = 1;
					break;
				}
			}
		}
		for (auto& c: suspcells){

			//edge intersection check
			for (auto& e: c->edges){
				if (e->id != 0) goto CONTINUE2;
			}
			ret.push_back(c);
CONTINUE2:
			continue;
		}
	}

	return ret;
}

CellData Grid::TExtractCells::extract_bound(const CellData& suspcells, const EdgeData& nted){
	CellData ret;
	EdgeData suspedges = AllEdges(suspcells);
	BoundingBox bb=HM2D::BBox(nted);

	aa::constant_ids_pvec(suspedges, 0);
	Finder::RasterizeEdges domainfinder(nted, bb, bb.maxlen()/30);
	double ksieta[2];
	for (auto& e: suspedges){
		Point p1 = *e->pfirst(), p2 = *e->plast();
		for (int isus: domainfinder.bbfinder().suspects(p1, p2)){
			Point *n1 = nted[isus]->pfirst(), *n2 = nted[isus]->plast();
			double m1 = Point::signed_meas_line(p1, *n1, *n2);
			double m2 = Point::signed_meas_line(p2, *n1, *n2);
			if (fabs(m1)<geps*geps){
				if (isOnSection(p1, *n1, *n2, ksieta[0])) {e->id = 1; break; }
			}
			if (fabs(m2)<geps*geps){
				if (isOnSection(p2, *n1, *n2, ksieta[0])) {e->id = 1; break; }
			}
			if (m1 < 0 && m2 < 0) continue;
			if (m1 > 0 && m2 > 0) continue;
			if (SectCross(*n1, *n2, p1, p2, ksieta)){ e->id = 1; break; }
		}
	}
	for (auto& c: suspcells){
		//edge intersection check
		for (auto& e: c->edges) if (e->id != 0){
			ret.push_back(c);
			break;
		}
	}

	return ret;
}

CellData Grid::TExtractCells::_run(const GridData& grid, const EdgeData& domain, int what){
	assert(Contour::IsContour(domain) && Contour::IsClosed(domain));
	Contour::Tree tree;
	tree.add_contour(domain);
	return ExtractCells(grid, tree, what);
}
