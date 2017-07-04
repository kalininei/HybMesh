#include "infogrid.hpp"
#include "contour_tree.hpp"
#include "modcont.hpp"
#include "finder2d.hpp"
#include "clipdomain.hpp"
#include "assemble2d.hpp"
#include "hmtimer.hpp"

using namespace HM2D;
using namespace HM2D::Finder;
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

class ExtractRunner{
	const int RasterN = 100;
	int what;
	CellData icells;
	EdgeData iedges;
	VertexData iverts;
	Contour::Tree ntdom;
	EdgeData dom_edges;
	vector<BoundingBox> gridbb;
	vector<BoundingBox> contbb;
	shared_ptr<RasterizeEdges> domraster;
	vector<int> raster_groups;
	vector<int> group_tp;   //0-undefined, 1-good, 2-bad
	vector<int> edge_pos;   //0-undefined, 1-good, 2-bad, 3-uncertain
	vector<int> cell_pos;   //0-undefined, 1-good, 2-bad, 3-uncertain, 4-fully uncertain
	
	bool sqrs_all_good(const vector<int>& sn){
		return std::all_of(sn.begin(), sn.end(),
				[&](int a){ return group_tp[raster_groups[a]] == 1; });
	}
	bool sqrs_all_bad(const vector<int>& sn){
		return std::all_of(sn.begin(), sn.end(),
				[&](int a){ return group_tp[raster_groups[a]] == 2; });
	}

	vector<vector<int>> _edge_squares;
	vector<bool> _has_edge_squares;
	vector<int>& edge_squares(int ied){
		if (!_has_edge_squares[ied]){
			_has_edge_squares[ied] = true;
			_edge_squares[ied] = domraster->bbfinder().sqrs_by_segment(
				*iedges[ied]->pfirst(), *iedges[ied]->plast());
		}
		return _edge_squares[ied];
	}

	CellData cells_by_tp(int tp){
		CellData ret;
		for (int i=0; i<icells.size(); ++i) if (cell_pos[i] == tp){
			ret.push_back(icells[i]);
		}
		return ret;
	}
public:
	CellData Result;

	ExtractRunner(const GridData& grid, const Contour::Tree& domain, int _what){
		what=_what;
		//get needed domain contours
		ntdom = domain;
		if (what != BOUND) ntdom.remove_detached();
		ntdom = Contour::Algos::Simplified(ntdom);
		//initial checks
		if (ntdom.nodes.size() == 0){
			if (what==OUTSIDE) Result=grid.vcells;
			return;
		}
		if (grid.vcells.size() == 0) return;

		//process bounding boxes and raterization.
		//fills all class internal data
		//#########
		icells = cells_by_tp(11);
		sort_by_bb(grid, domain);
		if (icells.size() == 0) return;

		//cell edges bounding box checks
		using_edge_bb();

		//cell edges cross checks
		using_edge_crosses();
	
		//cells which has intersections at edge nodes or
		//tangent intercections
		process_uncertain_cells();
	
		//check for domain areas which lie fully within cells
		detect_fully_inside_contours();

		//process cells with points liying on domain boundary (fully uncertain)
		if (what != BOUND) inner_point_check();
	
		//assemble result
		for (int i=0; i<icells.size(); ++i){
			if (cell_pos[i] != 2) Result.push_back(icells[i]);
		}
	}

	void sort_by_bb(const GridData& grid, const Contour::Tree& domain){
		//domain bounding box
		for (int i=0; i<ntdom.nodes.size(); ++i){
			contbb.push_back(HM2D::BBox(ntdom.nodes[i]->contour));
		}
		BoundingBox domainbb(contbb);

		//gridcells bounding box (+filtering cells out of domain bb)
		for (int i=0; i<grid.vcells.size(); ++i){
			auto& cell = grid.vcells[i];
			auto bb = HM2D::BBox(cell->edges);
			auto pos = domainbb.relation(bb);
			if (pos <= 1) {
				icells.push_back(cell);
				gridbb.push_back(bb);
			} else if (pos == 2 || pos == 4){
				//cell bb contain or intersects domain bb
				//this is either an outer or a bad cell
				if (what!=INSIDE){
					icells.push_back(cell);
					gridbb.push_back(bb);
				}
			} else if (pos == 3){
				//cell outside domain => good cell only for OUSIDE
				if (what!=INSIDE) Result.push_back(cell);
			}
		}
		if (icells.size() == 0) return;
		BoundingBox meshareabb(gridbb);

		//remove domain contours which lie outside mesharea domain
		std::set<int> to_remove;
		for (int i=0; i<ntdom.nodes.size(); ++i){
			auto rel = meshareabb.relation(contbb[i]);
			if (rel == 3) to_remove.insert(i);
		}
		assert(to_remove.size() < ntdom.nodes.size());
		if (to_remove.size() > 0){
			aa::remove_entries(ntdom.nodes, to_remove);
			aa::remove_entries(contbb, to_remove);
			domainbb = contbb[0];
			for (int i=1; i<contbb.size(); ++i){
				domainbb.widen(contbb[i]);
			}
		}

		//rasterization
		dom_edges = ntdom.alledges();
		domraster.reset(new RasterizeEdges(dom_edges, domainbb, domainbb.maxlen()/RasterN));

		raster_groups = domraster->colour_squares(what!=BOUND);
		if (what == BOUND){
			group_tp = std::vector<int>{0, 1};
		} else {
			int ngroups = *max_element(raster_groups.begin(),
						raster_groups.end()) + 1;
			group_tp.resize(ngroups, 0);
			for (int i=1; i<ngroups; ++i){
				int ns = std::find(raster_groups.begin(),
					raster_groups.end(), i)-raster_groups.begin();
				assert(ns < raster_groups.size());
				int pos = ntdom.whereis(
					domraster->bbfinder().sqr_center(ns));
				if (pos==what) group_tp[i]=1;
				else group_tp[i] = 2;
			}
		}

		//sort cells using rasterization
		for (int i=0; i<icells.size(); ++i){
			vector<int> sq = domraster->bbfinder().sqrs_by_bbox(
					gridbb[i]);
			assert(sq.size() != 0);
			//all are good
			if (sqrs_all_good(sq)){
				Result.push_back(icells[i]);
				icells[i].reset();
			}
			//all are bad
			if (sqrs_all_bad(sq)) icells[i].reset();
		}

		//restore all other data
		supplement_data();

		aa::constant_ids_pvec(grid.vcells, -1);
		aa::constant_ids_pvec(grid.vedges, -1);
		aa::constant_ids_pvec(grid.vvert, -1);
		aa::enumerate_ids_pvec(icells);
		aa::enumerate_ids_pvec(iedges);
		aa::enumerate_ids_pvec(iverts);
	}

	void supplement_data(){
		std::set<int> to_remove;
		for (int i=0; i<icells.size(); ++i)
			if (icells[i] == nullptr) to_remove.insert(i);
		aa::remove_entries(icells, to_remove);
		aa::remove_entries(gridbb, to_remove);
		iedges = AllEdges(icells);
		iverts=AllVertices(iedges);
		edge_pos.resize(iedges.size(), 0);
		cell_pos.resize(icells.size(), 0);
		_edge_squares.resize(iedges.size());
		_has_edge_squares.resize(iedges.size(), false);
		_vert_pos.resize(iverts.size(), -1);
		_dom_internal.resize(ntdom.nodes.size());
	}

	bool has_cross(const EdgeData& closed, const EdgeData& another){
		auto ca = Contour::Finder::CrossAll(another, closed);
		if (ca.size() == 0) return false;
		vector<double> candw; candw.reserve(ca.size()+2);
		for (int i=0; i<ca.size()-1; ++i){
			candw.push_back((std::get<2>(ca[i]) + std::get<2>(ca[i+1]))/2);
		}
		if (HM2D::Contour::IsClosed(another)){
			double a = (std::get<2>(ca[0]) + std::get<2>(ca.back()))/2.0;
			if (a>1) a-=1;
			candw.push_back(a);
		} else {
			candw.push_back( (std::get<2>(ca[0]))/2. );
			candw.push_back( (std::get<2>(ca.back()) + 1.0)/2. );
		}
		vector<Point> candp = Contour::WeightPoints(another, candw);

		for (int i=0; i<candp.size(); ++i){
			int pos = HM2D::Contour::Finder::WhereIs(closed, candp[i]);
			if (pos == INSIDE) return true;
		}
		return false;
	}
	
	//Sorting edge using bounding box and cross detections
	vector<int> _vert_pos;
	int vert_pos(int i){
		if (what == BOUND) return 1;
		if (_vert_pos[i] == -1){
			auto sq = domraster->bbfinder().sqrs_by_point(*iverts[i]);
			if (sq.size() == 0) _vert_pos[i] = what!=INSIDE ? 1 : 2;
			else if (sqrs_all_good(sq)) _vert_pos[i] = 1;
			else if (sqrs_all_bad(sq)) _vert_pos[i] = 2;
			else{
				int pos = ntdom.whereis(*iverts[i]);
				if (pos == BOUND) _vert_pos[i] = 0;
				else if (pos == what) _vert_pos[i] = 1;
				else _vert_pos[i] = 2;
			}
		}
		return _vert_pos[i];
	}

	vector<shared_ptr<Point>> _dom_internal;
	const Point& dom_internal(int ic){
		if (!_dom_internal[ic]){
			auto& cnt = ntdom.nodes[ic]->contour;
			Point cp;
			if (HM2D::Contour::IsClosed(cnt)){
				cp = Point(HM2D::Contour::InnerPoint(cnt));
			} else {
				cp = cnt[0]->center();
			}
			_dom_internal[ic].reset(new Point(cp));
		}
		return *_dom_internal[ic];
	}

	void good_edge(int i){
		edge_pos[i] = 1;
	};
	bool has_right_cell(int i){
		return iedges[i]->has_right_cell() && iedges[i]->right.lock()->id!=-1;
	}
	bool has_left_cell(int i){
		return iedges[i]->has_left_cell() && iedges[i]->left.lock()->id!=-1;
	}
	void uncertain_edge(int i){
		edge_pos[i] = 3;
		if (has_right_cell(i)) {
			if (cell_pos[iedges[i]->right.lock()->id] != 2)
				cell_pos[iedges[i]->right.lock()->id] = 3;
		}
		if (has_left_cell(i)) {
			if (cell_pos[iedges[i]->left.lock()->id] != 2)
				cell_pos[iedges[i]->left.lock()->id] = 3;
		}
	};
	void bad_edge(int i){
		edge_pos[i] = 2;
		if (has_right_cell(i)) {
			cell_pos[iedges[i]->right.lock()->id] = 2;
		}
		if (has_left_cell(i)) {
			cell_pos[iedges[i]->left.lock()->id] = 2;
		}
	};

	void analyze_edge_bbox(int i){
		if (edge_pos[i] != 0) return;
		auto edge = iedges[i];
		//if any end point is bad => edge is bad
		if (vert_pos(edge->pfirst()->id) == 2) return bad_edge(i);
		if (vert_pos(edge->plast()->id) == 2) return bad_edge(i);
		//if all finder squares are good => edge is good
		vector<int>& sqrs = edge_squares(i);
		if (sqrs_all_good(sqrs)) return good_edge(i);
		if (sqrs_all_bad(sqrs)) return bad_edge(i);
	};

	void analyze_edge_cross(int i){
		if (edge_pos[i] != 0) return;
		auto& ed = iedges[i];
		bool may_be_uncertain = false;
		for (int is: domraster->bbfinder().sqr_entries(edge_squares(i))){
			auto& ded = dom_edges[is];
			auto cr = SectCrossGeps(*ded->pfirst(), *ded->plast(), *ed->pfirst(), *ed->plast());
			if (cr.inner_cross()) return bad_edge(i);
			if (cr.has_contact()) may_be_uncertain = true;
		}
		if (may_be_uncertain) return uncertain_edge(i);
		else return good_edge(i);
	};

	bool domain_node_inside_cell(int icell, int icont, bool can_be_bound){
		const Point& wh = (can_be_bound) ? dom_internal(icont)
		                                 : *HM2D::Contour::First(ntdom.nodes[icont]->contour);
		auto pos = HM2D::Contour::Finder::WhereIs(icells[icell]->edges, wh);
		return pos == INSIDE;
	}

	void using_edge_bb(){
		//bad cells detection
		vector<bool> usede(iedges.size(), false);
		for (int i=0; i<icells.size(); ++i){
			auto cell = icells[i];
			int j=0;
			while (cell_pos[i] == 0 && j<cell->edges.size()){
				int ie = cell->edges[j++]->id;
				if (usede[ie]) continue;
				else usede[ie] = true;
				analyze_edge_bbox(ie);
			}
		}
	}

	void using_edge_crosses(){
		//bad and uncertain cells detection
		vector<bool> usede(iedges.size(), false);
		for (int i=0; i<icells.size(); ++i){
			auto cell = icells[i];
			int j=0;
			while (cell_pos[i] == 0 && j<cell->edges.size()){
				int ie = cell->edges[j++]->id;
				if (usede[ie]) continue;
				else usede[ie] = true;
				analyze_edge_cross(ie);
			}
		}
	}

	void detect_fully_inside_contours(){
		//undefined, uncertain => good or bad cells
		//Here all cells with pos = 0 have no edges which cross, are tangent to
		//domain contours or contain bad end points. So most likely these are
		//good cells except for cases when some domain contour lies fully inside those cells. 
		vector<int> last_level_dom;
		for (int i=0; i<ntdom.nodes.size(); ++i){
			if (ntdom.nodes[i]->children.size() == 0)
				last_level_dom.push_back(i);
		}
		for (int i=0; i<icells.size(); ++i) if (cell_pos[i] == 0 || cell_pos[i] >= 3){
			auto cell = icells[i];
			for (int j=0; j<last_level_dom.size(); ++j){
				int idom = last_level_dom[j];
				int bbpos = gridbb[i].relation(contbb[idom]);
				if (cell_pos[i]==4 && bbpos == 0){
					//equality:
					//not sure, but seems that above check is enough
					//since process_uncertain_cell wiped out all other cases.
					if ((ntdom.nodes[idom]->isinner() && what == INSIDE) ||
					    (ntdom.nodes[idom]->isouter() && what == OUTSIDE)){
						cell_pos[i] = 2;
					} else {
						cell_pos[i] = 1;
					}
					break;
				} else if (bbpos <= 1 && domain_node_inside_cell(i, idom, cell_pos[i]>=3)){
					cell_pos[i] = 2;
					break;
				}
			}
		}
	}

	void process_uncertain_cells(){
		//uncertain => uncertain, fully uncertain or bad
		//after that procedure cells marked as non-bad do not have
		//intersective crosses but only tangent ones.
		for (int i=0; i<icells.size(); ++i) if (cell_pos[i] == 3){
			//take all connected contours within cell box
			//and analyze position of points between possible crosses.
			//If any is inside then this is a bad cell.
			HM2D::EdgeData ae;
			for (int s: domraster->bbfinder().suspects(gridbb[i])){
				ae.push_back(dom_edges[s]);
			}
			for (auto& c: HM2D::Contour::Assembler::SimpleContours(ae)){
				if (has_cross(icells[i]->edges, c)){
					cell_pos[i] = 2;
					break;
				}
			}
		}
		//detect cells with all nodes lying on domain bnd.
		aa::enumerate_ids_pvec(iverts);
		if (what != BOUND)
		for (int i=0; i<icells.size(); ++i) if (cell_pos[i] == 3){
			cell_pos[i] = 4;
			for (auto& p: HM2D::Contour::OrderedPoints1(icells[i]->edges)){
				if (vert_pos(p->id) != 0){
					cell_pos[i] = 3;
					break;
				}
			}
		}
	}
	void inner_point_check(){
		for (int i=0; i<icells.size(); ++i) if (cell_pos[i] == 4){
			Point p = HM2D::Contour::InnerPoint(icells[i]->edges);
			int pos = ntdom.whereis(p);
			if ((pos == INSIDE && what == OUTSIDE) ||
			    (pos == OUTSIDE && what == INSIDE)){
				cell_pos[i] = 2;
			} else {
				cell_pos[i] = 1;
			}
		}
	}
};

}

HMCallback::FunctionWithCallback<Grid::TExtractCells> Grid::ExtractCells;

CellData Grid::TExtractCells::_run(const GridData& grid, const Contour::Tree& domain, int what){
	auto r = ExtractRunner(grid, domain, what);
	return r.Result;
}


CellData Grid::TExtractCells::_run(const GridData& grid, const EdgeData& domain, int what){
	assert(Contour::IsContour(domain) && Contour::IsClosed(domain));
	Contour::Tree tree;
	tree.add_contour(domain);
	return ExtractCells(grid, tree, what);
}
