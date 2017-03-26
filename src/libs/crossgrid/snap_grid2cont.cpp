#include "snap_grid2cont.hpp"
#include "finder2d.hpp"
#include "treverter2d.hpp"
#include "assemble2d.hpp"
#include "modgrid.hpp"
#include "healgrid.hpp"

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
			const VertexData& snap_nodes, bool only_corner_points=true): g(&grid){
		auto cv = HM2D::AllVertices(cont);
		//snapping nodes
		for (auto p: snap_nodes){
			//try to search among vertices
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
		//all contour significant points weights
		if (only_corner_points){
			cp = HM2D::Contour::CornerPoints1(cont);
		} else {
			cp = HM2D::Contour::OrderedPoints1(cont);
		}
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
		const VertexData& snap_nodes, bool only_significant){
	assert(Contour::IsContour(cont));
	//LeftCells reverter has index based algorithm. All new edges will be
	//added to the end of grid.vedges so it is save to used grid.vedges as
	//input data to the reverter.
	Contour::R::LeftCells grev(grid.vedges);
	Contour::R::ReallyDirect contrev(cont);

	auto proc = _ShiftSnapPreCalc(grid, cont, snap_nodes, only_significant);
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

GridData Algos::SnapBndSection(const GridData& grid, const EdgeData& cont,
		Vertex* gp1, Vertex* gp2, std::string snap_strategy){
	//checks
	assert(HM2D::Contour::IsContour(cont));
	bool isop = HM2D::Contour::IsOpen(cont);
	if ((gp1 == gp2 && isop) || (gp1 != gp2 && !isop)){
		throw std::runtime_error("open/close feature "
			"for grid boundary and contour should match");
	}
	//force clockwise direction of cont if it is closed
	HM2D::Contour::R::Clockwise cw(cont, false);

	//make a copy and reassign gp1, gp2
	GridData ret; HM2D::DeepCopy(grid, ret);
	aa::enumerate_ids_pvec(grid.vvert);
	gp1 = ret.vvert[gp1->id].get();
	gp2 = ret.vvert[gp2->id].get();

	//extract grid subcontour. Make a deep copy to safely revert
	auto gbnd = ECol::Assembler::GridBoundary(ret);
	HM2D::EdgeData subc;
	HM2D::DeepCopy(Contour::Assembler::Contour1(gbnd, gp1), subc, 0);
	Contour::R::Clockwise::Permanent(subc, false);
	if (gp1 != gp2){
		subc = Contour::Assembler::ShrinkContour(subc, gp1, gp2);
		Contour::R::ForceFirst::Permanent(subc, *Contour::First(cont));
	}

	//snap vertices
	vector<double> ww = Contour::EWeights(subc);
	vector<Point> np = Contour::WeightPoints(cont, ww);
	if (snap_strategy == "shift"){
		auto av = HM2D::AllVertices(cont);
		for (auto& p: np){
			auto fnd = HM2D::Finder::ClosestPoint(av, p);
			p.set(*av[std::get<0>(fnd)]);
		}
	}
	int i=0;
	for (auto p: Contour::OrderedPoints(subc)){
		p->set(np[i++]);
	}

	//add intermediate points.
	//revert points to guarantee (grid-is-left-to-contour) requirement.
	shared_ptr<Contour::R::ReallyRevert> rr;
	if (!subc[0]->has_left_cell()){
		rr.reset(new Contour::R::ReallyRevert(cont));
	}
	Algos::SnapToContour(ret, cont, {}, false);

	//check involved cells for counterclockwise rotation
	//do this before Heal cause heal will fix that.
	for (auto e: subc){
		HM2D::Cell* c = (e->has_left_cell()) ? e->left.lock().get()
						     : e->right.lock().get();
		if (Contour::Area(c->edges) < 0) 
			throw std::runtime_error("Resulting grid is invalid");
	}

	//shift may result in doubled vertices.
	//do heal before self intersection cross cause doubled points
	//will ruin SelfCross algorithm.
	if (snap_strategy == "shift") HM2D::Grid::Algos::Heal(ret);

	//check involved cells for self intersections
	for (auto e: subc){
		HM2D::Cell* c = (e->has_left_cell()) ? e->left.lock().get()
						     : e->right.lock().get();
		if (std::get<0>(Contour::Finder::SelfCross(c->edges)))
			throw std::runtime_error("Resulting grid is invalid");
	}

	return ret;
}
