#include "bgrid_impose.hpp"
#include "trigrid.h"
#include "contclipping.hpp"
#include "constructor.hpp"
#include "cont_assembler.hpp"
#include "algos.hpp"

using namespace HMBlay::Impl;

namespace{

class ContactArea{
	HM2D::Contour::Tree domain;
	const BGrid* grid;
	const int n_cfinder = 50; //partition for auxilliary grid
	vector<HM2D::EdgeData> intersections; //intersection contours
	std::map<int, int> intersect_results; //cell pair index -> intersection contour index
	std::set<int> cells_in_area; //all cells which were involved in intersections

	//unique integer fot two cells
	int __cellcell(const Cell* c1, const Cell* c2){
		int i1 = c1->get_ind(), i2 = c2->get_ind();
		if (i1>i2) std::swap(i1, i2);
		return i1*grid->n_cells()+i2;
	}

	void __set_cells_intersection(const Cell* c1, const Cell* c2){
		int ind = __cellcell(c1, c2);
		//if cells share the same feature they can not intersect
		if (grid->is_from_same_source(c1, c2)) return;
		//if result was already obtained
		if (intersect_results.find(ind) != intersect_results.end()) return;
		//apply intersection
		auto cont1 = Cell2Cont(c1);
		auto cont2 = Cell2Cont(c2);
		auto res = HM2D::Contour::Clip::Intersection(cont1, cont2);
		if (res.nodes.size() == 1){
			cells_in_area.insert(c1->get_ind());
			cells_in_area.insert(c2->get_ind());
			intersect_results[ind] = intersections.size();
			intersections.push_back(res.nodes[0]->contour);
		} else {
			intersect_results[ind] = -1;
		}
	}

	void FillDomain(){
		GGeom::Info::CellFinder cfinder(grid, n_cfinder, n_cfinder);
		for (int i=0; i<n_cfinder*n_cfinder; ++i){
			auto clist = cfinder.CellsBySquare(i);
			//check all pairs of cells positioned close to each other
			for (int j=0; j<clist.size(); ++j){
				for (int k=j+1; k<clist.size(); ++k){
					__set_cells_intersection(clist[j], clist[k]);
				}
			}
		}
		//unite all intersections
		domain = HM2D::Contour::Clip::Union(intersections);
		HM2D::Contour::Clip::Heal(domain);
	}

public:
	ContactArea(const BGrid* _grid): grid(_grid){
		FillDomain();
	}

	bool HasContact() const { return domain.nodes.size() > 0; }

	std::list<const Cell*> involved_cells() const {
		std::list<const Cell*> ret;
		for (auto i: cells_in_area) ret.push_back(grid->get_cell(i));
		return ret;
	}
};

HM2D::EdgeData
WidenCellCont(const HM2D::EdgeData& cont){
	vector<Point> ret;
	auto op = HM2D::Contour::CornerPoints1(cont);
	assert(op.size()>=4);
	if (op.size() < 4) return HM2D::EdgeData();
	vector<std::pair<Point, Point>> lines(cont.size());
	double h = 1e16;
	for (int i=0; i<(int)op.size()-1; ++i) {
		double h2 = Point::meas(*op[i], *op[i+1]);
		if (h2 < h) h = h2;
	}
	h = sqrt(h);
	for (int i=0; i<op.size()-1; ++i){
		Point* p1 = op[i].get();
		Point* p2 = op[i+1].get();
		double A[3] = { p1->y-p2->y, p2->x-p1->x, vecCrossZ(*p1,*p2) };
		Vect nrm {-A[0], -A[1]};
		vecSetLen(nrm, h);
		lines[i].first = *p1+nrm;
		lines[i].second = *p2+nrm;
	}
	double ksieta[2];
	for (int i=0; i<lines.size(); ++i){
		auto& ln  = lines[i];
		auto& ln1 = (i == lines.size()-1) ? lines[0] : lines[i+1];
		SectCrossWRenorm(ln.first, ln.second, ln1.first, ln1.second, ksieta);
		assert(ksieta[0]<gbig/2 && ksieta[1]<gbig/2);
		ret.push_back(Point::Weigh(ln.first, ln.second, ksieta[0]));
	}
	return HM2D::Contour::Constructor::FromPoints(ret, true);
}

bool CellsIntersect(const Cell* c1, const Cell* c2, const BGrid& grid1){
	//Finds if widened cells with different source feature are intersected 
	//Distance of widening is defined in WidenCellCont procedure.
	
	//if cells share the same feature return false
	if (grid1.is_from_same_source(c1, c2)) return false;

	//if have common edge return false
	for (int i=0; i<c1->dim(); ++i){
		int cp = 0;
		for (int j=0; j<c2->dim(); ++j){
			if (*c1->get_point(i) == *c2->get_point(j)) ++cp;
			if (cp>=2) return false;
		}
	}

	//calculate intersection of widened cells
	auto tcont1 = Cell2Cont(c1);
	auto tcont2 = Cell2Cont(c2);
	auto cont1 = WidenCellCont(tcont1);
	auto cont2 = WidenCellCont(tcont2);
	return std::get<0>(HM2D::Contour::Algos::Cross(cont1, cont2));
}

shared_ptr<BGrid> GridFromCells(shared_ptr<BGrid> grid, const std::list<const Cell*>& cells){
	auto ret = shared_ptr<BGrid>(new BGrid());
	if (cells.size() == 0) return ret;

	GGeom::Modify::ShallowAdd(grid.get(), ret.get());
	std::vector<bool> used(grid->n_cells(), false);
	for (auto c: cells) used[c->get_ind()] = true;
	std::vector<const Cell*> unused;
	for (int i=0; i<grid->n_cells(); ++i){
		const Cell* c = grid->get_cell(i);
		if (!used[c->get_ind()]) unused.push_back(c);
	}
	GGeom::Modify::RemoveCells(*ret, unused);
	return ret;
}

HM2D::Contour::Tree
AreaFromCells(const std::list<const Cell*>& plus, const std::list<const Cell*>& minus){
	vector<HM2D::EdgeData> cnts;
	//plus
	for (auto c: plus) cnts.push_back(Cell2Cont(c));
	auto r1 = HM2D::Contour::Clip::Union(cnts);
	//minus
	cnts.clear();
	for (auto c: minus) cnts.push_back(Cell2Cont(c));
	auto ret = HM2D::Contour::Clip::Difference(r1, cnts);
	HM2D::Contour::Clip::Heal(ret);
	return ret;
}


void TriAreaModify(HM2D::Contour::Tree& tree, const BGrid& og){
	//1) simplify
	tree = HM2D::Contour::Algos::Simplified(tree);
	//2) boundary points candidates
	auto bpv = GGeom::Info::BoundaryPoints(og);
	std::set<shared_ptr<GridPoint>> bp(bpv.begin(), bpv.end());
	BoundingBox bbox = HM2D::BBox(tree.alledges());
	auto it = bp.begin();
	while (it!=bp.end()){
		if (bbox.whereis(**it) == OUTSIDE) it = bp.erase(it);
		else ++it;
	}
	//3) place grid boundary points
	for (auto& n: tree.nodes){
		for (auto p: bp){
			auto fnd = HM2D::Contour::CoordAt(n->contour, *p);
			if (ISZERO(std::get<4>(fnd))){
				const double& w = std::get<3>(fnd);
				if (!ISZERO(w) && !ISEQ(w, 1)){
					HM2D::Contour::GuaranteePoint(n->contour, *p);
				}
			}
		}
	}
}

std::map<Point*, double> TriAreaWeights(const HM2D::Contour::Tree& tree){
	std::map<Point*, double> ret;
	for (auto& n: tree.nodes){
		vector<double> lens = HM2D::ELengths(n->contour);
		auto p = HM2D::Contour::OrderedPoints(n->contour);
		for (int i=0; i<p.size()-1; ++i){
			Point *pcur = p[i].get();
			double dist1 = (i==0) ? lens.back() : lens[i-1];
			double dist2 = lens[i];
			if (dist1>dist2) std::swap(dist1, dist2);
			//choose maximum not higher then coef*min
			ret[pcur] = std::min(dist2, 2.5*dist1);
		}
	}
	return ret;
}

//removes all grid area outside the closed cont
void PurgeGrid(BGrid& grid, const HM2D::EdgeData& cont){
	assert(HM2D::Contour::IsClosed(cont));
	bool innercont = HM2D::Contour::Area(cont) > 0;
	//fill bad points, good points: those which lie to the right of source contour
	//!! i'm still not sure whether tracking bad points to tell good/bad cells is enough.
	//Maybe there is a possibility of a not_too_bad_cell which is gathered all by bad points
	std::set<const GridPoint*> bad_points;
	for (int i=0; i<grid.n_points(); ++i){
		const GridPoint* p = grid.get_point(i);
		//Explicitly check points on contours because whereis is not relieable in this case
		auto res = HM2D::FindClosestEdge(cont, *p);
		if (ISZERO(std::get<1>(res))) continue;
		//check where is point
		int r = HM2D::Contour::WhereIs(cont, *p);
		if (innercont){
			if (r == OUTSIDE) bad_points.insert(p);
		} else {
			if (r == INSIDE) bad_points.insert(p);
		}
	}
	//fill bad cells: which contain bad_points (all points are bad/any point is bad)
	std::vector<const Cell*> bad_cells, not_too_bad_cells;
	for (int i=0; i<grid.n_cells(); ++i){
		const Cell* c = grid.get_cell(i);
		if (std::any_of(c->points.begin(), c->points.end(),
			[&](GridPoint* p){
				return bad_points.find(p) != bad_points.end();
			})){ bad_cells.push_back(c);}
		else continue;
		if (std::any_of(c->points.begin(), c->points.end(),
			[&](GridPoint* p){
				return bad_points.find(p) == bad_points.end();
			})){ not_too_bad_cells.push_back(c);}
	}

	//compensate good areas from not_too_bad_cells
	vector<HM2D::Contour::Tree> mgoodareas;
	vector<HM2D::EdgeData*> goodareas;
	for (auto c: not_too_bad_cells){
		auto ccont = Cell2Cont(c);
		HM2D::Contour::Clip::TRet ret;
		if (innercont) ret = HM2D::Contour::Clip::Intersection(ccont, cont);
		else ret = HM2D::Contour::Clip::Difference(ccont, cont);
		if (ret.nodes.size() < 1) continue;
		mgoodareas.push_back(ret);
		if (ret.nodes.size()>1) _THROW_NOT_IMP_
		else goodareas.push_back(&mgoodareas.back().nodes[0]->contour);
	}

	//remove bad cells
	GGeom::Modify::RemoveCells(grid, bad_cells);

	//compensate deleted areas
	int i=0;
	for (auto n: goodareas){
		vector<Point> p;
		for (auto pp: HM2D::Contour::CornerPoints(*n)) p.push_back(*pp);
		GGeom::Modify::AddCell(grid, p);
	}
}


//returns a closed contour from open one so that
//its closure edges not intersect grid area
const HM2D::EdgeData* ClosedSource(const HM2D::EdgeData* src, const BGrid* grid){
	if (HM2D::Contour::IsClosed(*src)) return src;
	//connecting first/last point of src so it doesn't cross grid area
	//basic contour for return. it is yet open.
	auto closedsrc = new HM2D::EdgeData();
	HM2D::DeepCopy(*src, *closedsrc);
	//assemble grid area
	auto sm = grid->subdivide();
	vector<HM2D::EdgeData> contvec;
	for (auto s: sm){
		auto tc = GGeom::Info::Contour(*s);
		for (auto r: tc.roots()) contvec.push_back(r->contour);
	}
	auto ar = HM2D::Contour::Clip::Union(contvec);
	HM2D::Contour::Clip::Heal(ar);

	//function to detect if section shortened by 4% intersects 
	//grid area. Shortage is needed to make end points negligible.
	auto goodline = [&ar, &closedsrc](Point& p1, Point& p2)->bool{
		auto w1 = std::make_shared<HM2D::Vertex>(Point::Weigh(p1, p2, 0.02));
		auto w2 = std::make_shared<HM2D::Vertex>(Point::Weigh(p1, p2, 0.98));
		HM2D::EdgeData wcont;
		wcont.emplace_back(new HM2D::Edge(w1, w2));
		//check cross with a source lint
		auto cres = HM2D::Contour::Algos::Cross(wcont, *closedsrc);
		if (std::get<0>(cres)) return false;
		//looping only over parent contours
		for (auto p: ar.roots()){
			auto cres = HM2D::Contour::Algos::Cross(wcont, p->contour);
			if (std::get<0>(cres)) return false;
		}
		return true;
	};
	//1) try: direct connection
	if (goodline(*HM2D::Contour::First(*src), *HM2D::Contour::Last(*src))){
		HM2D::Contour::AddLastPoint(*closedsrc, HM2D::Contour::First(*closedsrc));
		return closedsrc;
	}
	//2) try: connection via grid bounding box
	auto bbox1 = HM2D::BBox(ar.alledges()), bbox2 = HM2D::BBox(*src);
	BoundingBox bbox({bbox1, bbox2}); 
	bbox.widen(0.1*std::max(bbox.lenx(), bbox.leny()));
	auto sqrcont = HM2D::Contour::Constructor::FromPoints(
		{bbox.xmin, bbox.ymin, bbox.xmax, bbox.ymin,
		 bbox.xmax, bbox.ymax, bbox.xmin, bbox.ymax}, true);
	double sz = 1.1*sqrt(sqr(bbox.lenx()) + sqr(bbox.leny()));
	auto find_good_box_point = [&](Point& src)->shared_ptr<HM2D::Vertex>{
		static double an = 0;
		//trying different angles until we find non crossing section
		for (int i=0; i<20; ++i){
			Point p2 = src + Point(cos(an), sin(an)) * sz;
			if ((an += 1.0) > 2*M_PI) an -=2*M_PI;
			if (goodline(src, p2)){
				//find intersection of square and (src, p2) line
				auto linecont = HM2D::Contour::Constructor::FromPoints({src, p2});
				auto crossres = HM2D::Contour::Algos::Cross(sqrcont, linecont);
				return std::make_shared<HM2D::Vertex>(std::get<1>(crossres));
			}
		}
		return shared_ptr<HM2D::Vertex>();
	};
	auto epoint = find_good_box_point(*HM2D::Contour::Last(*closedsrc));
	if (epoint){
		//add epoint edge to source in order not to disregard possible crossing with newly added edge
		HM2D::Contour::AddLastPoint(*closedsrc, epoint);
		auto spoint = find_good_box_point(*HM2D::Contour::First(*closedsrc));
		if (spoint){
			//build subcontour from a square
			auto sqrsubcont = HM2D::Contour::Constructor::CutContour(sqrcont, *epoint, *spoint);
			//assemble and return
			auto p1 = HM2D::Contour::OrderedPoints(*closedsrc);
			auto p2 = HM2D::Contour::OrderedPoints(sqrsubcont);
			vector<Point> outp;
			for (auto p: p1) outp.push_back(*p);
			for (int i=1; i<p2.size(); ++i) outp.push_back(*p2[i]);
			delete closedsrc;
			return new HM2D::EdgeData(HM2D::Contour::Constructor::FromPoints(outp, true));
		}
	}
	//3) failed to close
	throw std::runtime_error("Can not find a way to close an open source contour");
}


void IntersectCellsInfo2(const BGrid& grid, std::list<const Cell*>& to_delete, std::list<const Cell*>& to_keep,
		std::function<double(const Cell*)> prifun){
	//build contact area
	ContactArea contact(&grid);
	if (!contact.HasContact()) return;

	//choose conflict cells
	std::list<const Cell*> confc = contact.involved_cells();
	
	//filter cells: cannot use ContactArea results
	//because they were obtained without widening
	auto cit = confc.begin();
	while (cit != confc.end()){
		//get all cells which intersect current
		vector<std::list<const Cell*>::iterator> adjcells;
		auto cit2 = confc.begin();
		while(cit2 != confc.end()){
			if (cit != cit2 && CellsIntersect(*cit, *cit2, grid)){
				adjcells.push_back(cit2);
			}
			++cit2;
		}
		//keep cit / keep cit2 / delete both
		double pri = prifun(*cit);
		bool keep_cit = true;
		for (auto i: adjcells){
			double pri2 = prifun(*i);
			if (ISEQGREATER(pri2, pri)){
				keep_cit = false; break;
			}
		}
		//adding to list
		if (!keep_cit) to_delete.push_back(*cit);
		else to_keep.push_back(*cit);
		++cit;
	}
}

shared_ptr<BGrid> NonIntersectingGrid(const BGrid& grid, std::list<const Cell*>& to_delete, std::list<const Cell*>& to_keep){
	//build grid with old cells
	shared_ptr<BGrid> ret(new BGrid);
	ret->ShallowAdd(grid);
	GGeom::Modify::RemoveCells(*ret, vector<const Cell*>(to_delete.begin(), to_delete.end()));

	//build an area from deleted cells
	auto delarea = AreaFromCells(to_delete, to_keep);
	
	//modify delarea: simplify + place boundary points 
	TriAreaModify(delarea, *ret);

	//calculates weights for triangulation
	std::map<Point*, double> triweights = TriAreaWeights(delarea);
	
	//triangulate
	auto fillg = TriGrid::TriangulateArea(delarea, triweights, 100);
	//remove zero size triangles which may appear due to singular areas in delarea
	auto ars3 = fillg->cell_areas();
	vector<const Cell*> rmcells;
	for (int i=0; i<ars3.size(); ++i) if (ars3[i]<geps*geps) rmcells.push_back(fillg->get_cell(i));
	GGeom::Modify::RemoveCells(*fillg, rmcells);

	//add triangle grid and return
	GGeom::Modify::ShallowAdd(fillg.get(), ret.get());

	//heal grid, get rid of hanging nodes
	GGeom::Repair::Heal(*ret);
	NoHangingNodes nhn(*ret);

	return ret;
}

}

void NoHangingNodes::FindCellForNode(GridPoint* p, Cell*& c, int& ind){
	auto candv = cfinder.CellCandidates(*p);
	double ksi;
	for (auto cand: candv){
		for (int i=0; i<cand->dim(); ++i){
			isOnSection(*p, *cand->get_point(i),
					*cand->get_point(i+1), ksi);
			if (ksi>geps && ksi<1-geps){
				ind = i;
				c = const_cast<Cell*>(cand);
				return;
			}
		}
	}
	//not found. Perhaps p was put into hanging nodes list by mistake
	c = 0;
	ind = -1;
}
void NoHangingNodes::PlaceNode(GridPoint* p, Cell* c, int ind){
	while (ind>=c->dim()) ind-=c->dim();
	ind = (ind == c->dim()-1) ? 0 : ind+1;
	c->points.insert(c->points.begin() + ind, p);
}

void NoHangingNodes::CellSplit(Cell* cell, const GridPoint* p1, const GridPoint* p2){
	auto ind1 = std::find(cell->points.begin(), cell->points.end(), p1)
		- cell->points.begin();
	auto ind2 = std::find(cell->points.begin(), cell->points.end(), p2)
		- cell->points.begin();
	//build cells
	if (ind1>ind2) std::swap(ind1, ind2);
	shared_ptr<Cell> cell1(new Cell), cell2(new Cell);
	for (int i=ind1; i<=ind2; ++i){
		GridPoint* gp = const_cast<GridPoint*>(cell->get_point(i));
		cell1->points.push_back(gp);
	}
	for (int i=ind2; i<=ind1 + cell->dim(); ++i){
		GridPoint* gp = const_cast<GridPoint*>(cell->get_point(i));
		cell2->points.push_back(gp);
	}
	cell->points = cell1->points;

	//add new cell
	grid->ShallowAddCell(cell2, cell);
	cfinder.AddCell(cell2.get());
}

void NoHangingNodes::SimpleSplit(Cell* cell, int ind){
	//find closest cell node which doesn't
	//lie on the same edge as ind node
	auto p0 = cell->get_point(ind-1);
	auto p1 = cell->get_point(ind);
	auto p2 = cell->get_point(ind+1);
	Vect nextv = *p2 - *p1; vecNormalize(nextv);
	Vect prevv = *p0 - *p1; vecNormalize(prevv);

	int ind2 = -1; double meas = 1e20;
	for (int i=ind+2; i<cell->dim()+ind-1; ++i){
		Vect v = *cell->get_point(i) - *p1; vecNormalize(v);
		if (v == prevv || v == nextv) continue;
		double m = Point::meas(*p1, *cell->get_point(i));
		if (ind2<0 || m<meas){
			ind2 = i;
			meas = m;
		}
	}
	while (ind2>=cell->dim()) ind2-=cell->dim();
	assert(ind2>=0);

	//resulting cells have lowest priority since they are no more regular
	grid->RemoveFeatures(cell);

	CellSplit(cell, p1, cell->get_point(ind2));
}

void NoHangingNodes::SimplifyCell(Cell* cell, int ind){
	while (ind>=cell->dim()) ind-=cell->dim();
	auto pind = cell->get_point(ind);
	if (cell->dim() == 4){
		auto pind2 = cell->get_point(ind+2);
		CellSplit(cell, pind, pind2);
	} else if (cell->dim() == 5 && grid->is_from_source(cell)){
		//find opposite node
		double w = Point::dist(*pind, *cell->get_point(ind-1))/
			Point::dist(*cell->get_point(ind+1), *cell->get_point(ind-1));
		Point tp = Point::Weigh(*cell->get_point(ind+2), *cell->get_point(ind+3), 1-w);
		shared_ptr<GridPoint> p(new GridPoint(tp));
		//find another cell containing this node
		const Cell* cell2 = cfinder.FindExcept(*p, {cell});
		//analyze cell2
		if (cell2 == 0){
			//if no such cell just add a node to boundary
			grid->ShallowAddNode(p);
			PlaceNode(p.get(), cell, ind+2);
			CellSplit(cell, pind, p.get());
		} else if (grid->get_weight(cell) - grid->get_weight(cell2) == 1){
			//if cell2 is on lower layer split it
			grid->ShallowAddNode(p);
			PlaceNode(p.get(), cell, ind+2);
			CellSplit(cell, pind, p.get());
			AddHangingNode(p.get());
		} else {
			//if layers are same use something different
			_DUMMY_FUN_;
			SimpleSplit(cell, ind);
		}
	} else SimpleSplit(cell, ind);
}

void NoHangingNodes::AddHangingNode(GridPoint* p){
	Cell* cell;
	int ind;
	FindCellForNode(p, cell, ind);
	if (cell){
		PlaceNode(p, cell, ind);
		SimplifyCell(cell, ind+1);
	}
}

void NoHangingNodes::AnalyzeSingularContour(
		const HM2D::EdgeData& c, std::list<GridPoint*>& lst){
	//put all points except corner points to lst
	auto info = HM2D::Contour::OrderedInfo(c);
	for (int i=0; i<c.size(); ++i){
		Point* ps[3] {info[i].pprev.get(), info[i].p.get(), info[i].pnext.get()};
		Vect v1 = *ps[0] - *ps[1]; vecNormalize(v1);
		Vect v2 = *ps[2] - *ps[1]; vecNormalize(v2);
		if (Point::meas(v1,v2)>geps){
			lst.push_back(static_cast<GridPoint*>(ps[1]));
		}
	}
}
std::list<GridPoint*> NoHangingNodes::Find(){
	std::list<GridPoint*> hnodes;
	//hanging nodes create zero-area domains in grid contour tree
	auto gcont = GGeom::Info::Contour(*grid);
	for (auto n: gcont.nodes){
		bool status = false;
		auto info = HM2D::Contour::OrderedInfo(n->contour);
		for (int i1 = 0; i1<n->contour.size(); ++i1){
			Point* ps[3] {info[i1].pprev.get(), info[i1].p.get(), info[i1].pnext.get()};
			double ar = triarea(*ps[0], *ps[1], *ps[2]);
			if (fabs(ar)>geps) {status = true; break;}
		}
		if (!status){
			AnalyzeSingularContour(n->contour, hnodes);
		}
	}
	return hnodes;
}

NoHangingNodes::NoHangingNodes(BGrid& _grid): grid(&_grid), cfinder(grid, 20, 20){
	//build list of hanging nodes
	auto hnodes = Find();
	//modify each
	for (auto p: hnodes) AddHangingNode(p);
}


shared_ptr<BGrid> HMBlay::Impl::BGridImpose(shared_ptr<BGrid> grid,
		std::function<double(const Cell*)> prifun,
		const HM2D::EdgeData& source){
	//1) if source is not a closed contour supplement it
	const HM2D::EdgeData* src = ClosedSource(&source, grid.get());

	//2) remove cells or their parts which lie to the right side of source
	PurgeGrid(*grid, *src);

	//don't need it anymore
	if (src != &source) delete src;

	//3) analyse grid for intersecting cells.
	//   returns two sets of intersecting cells: those which should be deleted and keeped
	std::list<const Cell*> to_delete, to_keep;
	IntersectCellsInfo2(*grid, to_delete, to_keep, prifun);
	if (to_delete.size() == 0) return grid;
	
	//4) build a grid: non-intersecting cells + to_keep_cells + triangled to_delete area
	return NonIntersectingGrid(*grid, to_delete, to_keep);
}

HM2D::EdgeData HMBlay::Impl::Cell2Cont(const Cell* c){
	HM2D::VertexData cp;
	for (auto i=0; i<c->dim(); ++i){
		cp.push_back(std::make_shared<HM2D::Vertex>(*c->get_point(i)));
	}
	return HM2D::Contour::Assembler::Contour1(cp, true);
}
