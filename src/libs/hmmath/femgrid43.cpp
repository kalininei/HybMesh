#include "femgrid43.hpp"
#include "trigrid.h"
#include "procgrid.h"
#include "debug_grid2d.h"

using namespace HMFem;

bool Grid43::check(){
	GGeom::Debug::save_vtk(*this, "before.vtk");
	for (int i=0; i<n_cells(); ++i){
		if (cells[i]->dim() == 3){
			if (triarea(*cells[i]->points[0],
					*cells[i]->points[1],
					*cells[i]->points[2]) < -geps*geps)
				return false;
		} else {
			_THROW_NOT_IMP_;
		}
	}
	return true;
}

shared_ptr<Grid43>
Grid43::Build(GridGeom* parent){
	shared_ptr<Grid43> ret(new Grid43());
	//Temp. implementation: only 4/3 grids on input
	assert([](GridGeom* g){
		for (int i=0; i<g->n_cells(); ++i){
			if (g->get_cell(i)->dim()>4) return false;
		}
		return true;
	}(parent));
	//shallow copy of points, deep copy of cells
	shallow_copy(parent, ret.get());

	ret->cells.clear();
	for (int i=0; i<parent->n_cells(); ++i){
		auto cell_old = parent->get_cell(i);
		auto c = aa::add_shared(ret->cells, Cell());
		for (int j=0; j<cell_old->dim(); ++j){
			add_point_to_cell(c, const_cast<GridPoint*>(cell_old->get_point(i)));
		}
	}

	return ret;
}

shared_ptr<Grid43>
Grid43::Build3(const vector<Point>& p, double h){
	shared_ptr<Grid43> ret(new Grid43());
	shared_ptr<TriGrid> g1 = TriGrid::TriangulateArea(p, h);
	GGeom::Modify::ShallowAdd(g1.get(), ret.get());
	return ret;
}

shared_ptr<Grid43>
Grid43::Build3(const vector<vector<Point>>& cnts, double h){
	shared_ptr<Grid43> ret(new Grid43());
	shared_ptr<TriGrid> g1 = TriGrid::TriangulateArea(cnts, h);
	GGeom::Modify::ShallowAdd(g1.get(), ret.get());
	return ret;
}

shared_ptr<Grid43>
Grid43::Build3(const HMCont2D::ContourTree& conts,
		const ShpVector<HMCont2D::Contour>& constr,
		double h){
	TriGrid g1(conts, constr, h);
	shared_ptr<Grid43> ret(new Grid43());
	GGeom::Modify::ShallowAdd(&g1, ret.get());
	return ret;
}
shared_ptr<Grid43>
Grid43::Build3(const HMCont2D::ContourTree& conts,
	const ShpVector<HMCont2D::Contour>& constraints,
	std::map<Point*, double>& h, double hh){
	TriGrid g1(conts, constraints, h, hh);
	shared_ptr<Grid43> ret(new Grid43());
	GGeom::Modify::ShallowAdd(&g1, ret.get());
	return ret;
}

shared_ptr<Grid43>
Grid43::Build3(GridGeom* orig){
	shared_ptr<Grid43> ret(new Grid43());
	//copy points and cells
	GGeom::Modify::ShallowAdd(orig, ret.get());
	//remove cells
	ret->cells.clear();
	//build triangle cells
	auto addtri = [&](Cell* oldcell, int i1, int i2, int i3){
		Cell* newcell = aa::add_shared(ret->cells, Cell());
		add_point_to_cell(newcell, const_cast<GridPoint*>(oldcell->get_point(i1)));
		add_point_to_cell(newcell, const_cast<GridPoint*>(oldcell->get_point(i2)));
		add_point_to_cell(newcell, const_cast<GridPoint*>(oldcell->get_point(i3)));
	};
	for (int i=0; i<orig->n_cells(); ++i){
		auto oc = orig->get_cell(i);
		if (oc->dim() == 3) addtri(oc, 0, 1, 2);
		else if (oc->dim() == 4){ addtri(oc, 0, 1, 2); addtri(oc, 0, 2, 3); }
		else _THROW_NOT_IMP_;
	}
	//set indicies
	ret->set_indicies();
	return ret;
}

void Grid43::AddSegments(const vector<vector<Point>>& pts){
	if (pts.size()==0) return;
	auto eds = get_edges();
	//find cells
	std::map<int, int> cells1;  //cells with one edge to divide
	std::map<int, vector<int>> cells2; //cells with multiple edges to divide
	ShpVector<GridPoint> newpoints;
	vector<vector<GridPoint*>> edpnt(pts.size());
	auto add_edge_cell = [&](int i, int icell){
		auto mapfnd = cells1.find(icell);
		if (mapfnd == cells1.end()) cells1[icell] = i;
		else {
			if (cells2.find(icell) == cells2.end()){
				cells2.emplace(icell, vector<int> {mapfnd->second, i});
				cells1.erase(mapfnd);
			} else cells2[icell].push_back(i);
		}
	};
	for (int i=0; i<pts.size(); ++i){
		auto& pe = pts[i];
		auto fnd1 = std::find_if(points.begin(), points.end(),
				[&pe](shared_ptr<GridPoint>& p){ return *p == pe[0]; });
		auto fnd2 = std::find_if(points.begin(), points.end(),
				[&pe](shared_ptr<GridPoint>& p){ return *p == pe.back(); });
		assert(fnd1 != points.end() && fnd2 != points.end());
		int p1 = fnd1 - points.begin();
		int p2 = fnd2 - points.begin();
		//new points
		edpnt[i].push_back(points[p1].get());
		for (int j=1; j<pts[i].size()-1; ++j){
			newpoints.emplace_back(new GridPoint(pts[i][j]));
			edpnt[i].push_back(newpoints.back().get());
		}
		edpnt[i].push_back(points[p2].get());
		//edges
		auto efnd = eds.find(Edge(p1, p2));
		assert(efnd != eds.end());
		if (efnd->cell_left>=0) add_edge_cell(i, efnd->cell_left);
		if (efnd->cell_right>=0) add_edge_cell(i, efnd->cell_right);
	}
	//construct new cells
	ShpVector<Cell> newcells;
	//treat cells with single boundary edge
	for (auto& it: cells1){
		auto& pe = pts[it.second];
		auto& pp = edpnt[it.second];
		auto& c = cells[it.first];
		//find first point
		int p0 = 0;
		while (c->points[p0] != pp[0]) ++p0;
		bool isleft = (c->get_point(p0+1) == pp.back());
		if (!isleft) { std::reverse(pp.begin(), pp.end()); p0 -= 1;}
		//create new cells
		for (int j=0; j<pp.size()-1; ++j){
			Cell* nc = aa::add_shared(newcells, Cell());
			nc->points.resize(3);
			nc->points[0] = const_cast<GridPoint*>(c->get_point(p0-1));
			nc->points[1] = pp[j];
			nc->points[2] = pp[j+1];
		}
	}
	//cells with multiple boundary edges
	for (auto& it: cells2){
		Cell& c = *cells[it.first];
		vector<GridPoint*> ppts;
		for (int i=0; i<c.dim(); ++i){
			GridPoint* p1 = c.points[i];
			GridPoint* p2 = c.points[(i+1) % c.dim()];
			ppts.push_back(p1);
			for (int k=0; k<it.second.size(); ++k){
				auto& pp = edpnt[it.second[k]];
				GridPoint* start = pp[0];
				GridPoint* end = pp.back();
				if (p1 == start && p2 == end){
					ppts.insert(ppts.end(), pp.begin()+1, pp.end()-1);
					break;
				} else if (p1 == end && p2 == start){
					ppts.insert(ppts.end(), pp.rbegin()+1, pp.rend()-1);
					break;
				}
			}
		}
		double h=0;
		for (int i=0; i<ppts.size(); ++i){
			double hh = Point::meas(*ppts[i], *ppts[(i+1) % ppts.size()]);
			if (hh > h) h = hh;
		}
		vector<Point> pts(ppts.size());
		for (int i=0; i<ppts.size(); ++i) pts[i] = *ppts[i];
		auto trigrid = TriGrid::TriangulateArea(pts, sqrt(h));
		for (int i=0; i<trigrid->n_cells(); ++i){
			Cell* c3 = trigrid->get_cell(i);
			Cell* nc = aa::add_shared(newcells, Cell());
			for (int j=0; j<c3->dim(); ++j){
				GridPoint  p = *c3->get_point(j);
				auto fnd = std::find_if(ppts.begin(), ppts.end(),
						[&p](GridPoint* x){ return *x == p; });
				if (fnd == ppts.end()){
					ppts.push_back(aa::add_shared(newpoints, p));
					nc->points.push_back(ppts.back());
				} else nc->points.push_back(*fnd);
			}
		}
	}
	//reconstruct grid
	std::set<int> cells_to_remove;
	for (auto& it: cells1) cells_to_remove.insert(it.first);
	for (auto& it: cells2) cells_to_remove.insert(it.first);
	aa::remove_entries(cells, cells_to_remove);
	cells.insert(cells.end(), newcells.begin(), newcells.end());
	points.insert(points.end(), newpoints.begin(), newpoints.end());
	set_indicies();
}

shared_ptr<Grid43::Approximator> Grid43::GetApprox() const{
	return shared_ptr<Grid43::Approximator>(new Grid43::Approximator(this));
}

Point Grid43::Approximator::LocalCoordinates(const Cell* c, Point p){
	if (c->dim() == 3) return LocalCoordinates3(c, p);
	else return LocalCoordinates4(c, p);
}

Point Grid43::Approximator::LocalCoordinates3(const Cell* c, Point p){
	auto &x1 = c->get_point(0)->x, &x2 = c->get_point(1)->x, &x3 = c->get_point(2)->x;
	auto &y1 = c->get_point(0)->y, &y2 = c->get_point(1)->y, &y3 = c->get_point(2)->y;

	double j11 = x2 - x1, j21 = x3 - x1;
	double j12 = y2 - y1, j22 = y3 - y1;
	double modj = (j22*j11 - j21*j12);

	if (fabs(modj) < geps*geps){
		double ksi;
		if (isOnSection(p, *c->get_point(0), *c->get_point(1), ksi)){
			return Point(ksi, 0);
		} else if (isOnSection(p, *c->get_point(1), *c->get_point(2), ksi)){
			return Point(1-ksi, ksi);
		} else if (isOnSection(p, *c->get_point(0), *c->get_point(2), ksi)){
			return Point(0, ksi);
		} else {
			assert(false);
		}
	}

	Point ksieta;
	ksieta.x = ( j22*(p.x - x1) - j21*(p.y - y1))/modj;
	ksieta.y = (-j12*(p.x - x1) + j11*(p.y - y1))/modj;

	//
	//assert( ksieta.x > -0.1 && ksieta.x < 1.1 );
	//assert( ksieta.y > -0.1 && ksieta.y < 1-ksieta.x+0.1 );
	//

	return ksieta;
}

Point Grid43::Approximator::LocalCoordinates4(const Cell* c, Point p){
	_THROW_NOT_IMP_;
}

double Grid43::Approximator::Interpolate(const Cell* c, Point ksieta, const vector<double>& fun){
	if (c->dim() == 3) return Interpolate3(c, ksieta, fun);
	else return Interpolate4(c, ksieta, fun);
}
double Grid43::Approximator::Interpolate3(const Cell* c, Point ksieta, const vector<double>& fun){
	int i1 = c->get_point(0)->get_ind(),
	    i2 = c->get_point(1)->get_ind(),
	    i3 = c->get_point(2)->get_ind();
	return (1-ksieta.x-ksieta.y)*fun[i1] + (ksieta.x)*fun[i2] + (ksieta.y)*fun[i3];
}

double Grid43::Approximator::Interpolate4(const Cell* c, Point ksieta, const vector<double>& fun){
	_THROW_NOT_IMP_;
}

void Grid43::Approximator::FillJ3(std::array<double, 5>& J, const Cell* c){
	auto &x1 = c->get_point(0)->x, &x2 = c->get_point(1)->x, &x3 = c->get_point(2)->x;
	auto &y1 = c->get_point(0)->y, &y2 = c->get_point(1)->y, &y3 = c->get_point(2)->y;

	J[1] = x2 - x1; J[3] = x3 - x1;
	J[2] = y2 - y1, J[4] = y3 - y1;
	J[0] = (J[4]*J[1] - J[3]*J[2]);
}

void Grid43::Approximator::FillJ4(std::array<double, 5>& J, const Point& p, const Cell* c){
	_THROW_NOT_IMP_;
}

const Cell* Grid43::Approximator::FindPositive(const Point& p, Point& ksieta) const{
	auto candidates = cfinder->CellCandidates(p);
	std::array<double, 5> J; //modj, j11, j12, j21, j22
	std::vector<Point> bad_ksieta(candidates.size());
	for (int i=0; i<candidates.size(); ++i){
		auto& c = candidates[i];
		auto& ke = bad_ksieta[i];
		if (c->dim() == 3){FillJ3(J, c);}
		else if (c->dim() == 4){ _THROW_NOT_IMP_;}
		else { assert(false); }
		if (fabs(J[0]) < geps*geps){
			double ksi;
			if (isOnSection(p, *c->get_point(0), *c->get_point(1), ksi)){
				ke = Point(ksi, 0);
			} else if (isOnSection(p, *c->get_point(1), *c->get_point(2), ksi)){
				ke = Point(1-ksi, ksi);
			} else if (isOnSection(p, *c->get_point(0), *c->get_point(2), ksi)){
				ke = Point(0, ksi);
			} else assert(false);
		} else {
			auto cp = c->get_point(0);
			ke.x = ( J[4]*(p.x - cp->x) - J[3]*(p.y - cp->y))/J[0];
			ke.y = (-J[2]*(p.x - cp->x) + J[1]*(p.y - cp->y))/J[0];
		}
		if (J[0] > -geps*geps && ke.x>-geps && ke.x<1+geps && ke.y>-geps && ke.y<1-ke.x+geps){
			ksieta = ke;
			return c;
		}
	}
	//find best outer approximation
	int besti=-1;
	double bestdist = std::numeric_limits<double>::max();
	for (int i=0; i<candidates.size(); ++i){
		auto& ke = bad_ksieta[i];
		//if negative triangle contains point
		if (ke.x>-geps && ke.x<1+geps && ke.y>-geps && ke.y<1-ke.x+geps){
			besti = i;
			break;
		}
		//point lies outside a positive triangle: find measure to the basis triangle
		double d = std::min(Point::meas_section(ke, Point(0,0), Point(1,0)),
			std::min(Point::meas_section(ke, Point(1,0), Point(0,1)),
				Point::meas_section(ke, Point(0,1), Point(0,0))));
		if (d<bestdist) besti = i;
	}
	if (besti>-1){
		ksieta = bad_ksieta[besti];
		return candidates[besti];
	} else throw GGeom::EOutOfArea(p);
}

double Grid43::Approximator::Val(Point p, const vector<double>& fun) const{
	Point ksieta;
	const Cell* c = FindPositive(p, ksieta);
	return Interpolate(c, ksieta, fun);
}


vector<double> Grid43::Approximator::Vals(Point p, const vector<const vector<double>*>& funs) const{
	Point ksieta;
	const Cell* c = FindPositive(p, ksieta);
	vector<double> ret;
	for (auto& fun: funs) ret.push_back(Interpolate(c, ksieta, *fun));
	return ret;
}

const vector<Edge>& Grid43::Approximator::BndEdgesByPnt(const Point& p) const{
	if (bndedges.size() == 0){
		bndedges.resize(cfinder->NumSquares(), vector<Edge>());
		std::set<Edge> ed = grid->get_edges();
		for (auto e: ed) if (e.is_boundary()){
			BoundingBox bb(*grid->get_point(e.p1), *grid->get_point(e.p2));
			BoundingBox bb2(bb);
			bb2.xmin -= 0.1*(bb.xmax - bb.xmin);
			bb2.xmax += 0.1*(bb.xmax - bb.xmin);
			bb2.ymin -= 0.1*(bb.ymax - bb.ymin);
			bb2.ymax += 0.1*(bb.ymax - bb.ymin);
			auto ind = cfinder->IndSet(bb2);
			for (auto i: ind) if (i>=0 && i<cfinder->NumSquares()){
				bndedges[i].push_back(e);
			}
		}
	}
	return bndedges[cfinder->GetSquare(p)];
}
double Grid43::Approximator::BndVal(const Point& p, const vector<double>& fun) const{
	vector<Edge> ec = BndEdgesByPnt(p);
	vector<double> ksi(ec.size());
	vector<double> dist(ec.size());
	for (int i=0; i<ksi.size(); ++i){
		auto& edge = ec[i];
		dist[i] = Point::meas_section(p, *grid->get_point(edge.p1), *grid->get_point(edge.p2), ksi[i]);
	}
	int mi = std::min_element(dist.begin(), dist.end()) - dist.begin();
	assert(mi<dist.size());
	auto& edge = ec[mi];
	return (1 - ksi[mi])*fun[edge.p1] + ksi[mi]*fun[edge.p2];
}
vector<double> Grid43::Approximator::BndVals(const Point& p, const vector<const vector<double>*>& funs) const{
	vector<Edge> ec = BndEdgesByPnt(p);
	vector<double> ksi(ec.size());
	vector<double> dist(ec.size());
	for (int i=0; i<ksi.size(); ++i){
		auto& edge = ec[i];
		dist[i] = Point::meas_section(p, *grid->get_point(edge.p1), *grid->get_point(edge.p2), ksi[i]);
	}
	int mi = std::min_element(dist.begin(), dist.end()) - dist.begin();
	assert(mi<dist.size());
	auto& edge = ec[mi];
	vector<double> ret;
	for (auto& fun: funs){
		double v = (1 - ksi[mi])*(*fun)[edge.p1] + ksi[mi]*(*fun)[edge.p2];
		ret.push_back(v);
	}
	return ret;
}

// ============================ Isoline Builder
Grid43::IsolineBuilder::IsolineBuilder(shared_ptr<Grid43>& grid43,
		shared_ptr<Grid43::Approximator> approx43){
	//force grid to have only triangle cells
	bool only3 = true;
	for (int i=0; i<grid43->n_cells(); ++i)
		if (grid43->get_cell(i)->dim() == 4){
			only3 = false; break;
		}
	if (only3){
		grid = grid43;
		if (approx43) approx = approx43;
		else approx = grid->GetApprox();
	} else {
		_THROW_NOT_IMP_;
	}
}

bool Grid43::IsolineBuilder::AddLine(const Cell* c,
		HMCont2D::ECollection& ecol,
		HMCont2D::PCollection& pcol,
		double value, const vector<double>& fun) const{
	int i0 = c->get_point(0)->get_ind();
	int i1 = c->get_point(1)->get_ind();
	int i2 = c->get_point(2)->get_ind();
	double f0 = fun[i0], f1 = fun[i1], f2 = fun[i2];
	if (ISEQ(f0, f1) && ISEQ(f1, f2)) return false;
	//sort f0<f1<f2
	if (f1>f2) {std::swap(i1, i2); std::swap(f1, f2); };
	if (f0>f1) {std::swap(i0, i1); std::swap(f0, f1); };
	if (f1>f2) {std::swap(i1, i2); std::swap(f1, f2); };

	//check all possible positions
	if (value<f0 || value > f2) return false;
	else if (ISEQ(value, f0)){
		if (ISGREATER(f1, f0)) return false;
		else BuildLine(i0, i1, 0, i0, i1, 1, ecol, pcol);
	} else if (ISEQ(value, f1)){
		if (!ISEQ(f1,f2)) return false;
		else BuildLine(i1, i2, 0, i1, i2, 1, ecol, pcol);
	} else if (ISEQ(value, f2)){
		return false;
	} else if (value < f1){
		double w1 = (value - f0)/(f1 - f0);
		double w2 = (value - f0)/(f2 - f0);
		BuildLine(i0, i1, w1, i0, i2, w2, ecol, pcol);
	} else if (value < f2){
		double w1 = (value - f1)/(f2 - f1);
		double w2 = (value - f0)/(f2 - f0);
		BuildLine(i1, i2, w1, i0, i2, w2, ecol, pcol);
	} else return false;

	return true;
}

void Grid43::IsolineBuilder::BuildLine(int i0, int i1, double wi,
		int j0, int j1, double wj,
		HMCont2D::ECollection& ecol,
		HMCont2D::PCollection& pcol) const{
	auto pi = std::make_shared<Point>(
		Point::Weigh(*grid->get_point(i0), *grid->get_point(i1), wi));
	auto pj = std::make_shared<Point>(
		Point::Weigh(*grid->get_point(j0), *grid->get_point(j1), wj));
	pcol.add_value(pi);
	pcol.add_value(pj);
	ecol.add_value(HMCont2D::Edge(pi.get(), pj.get()));
}

HMCont2D::Container<HMCont2D::Contour>
Grid43::IsolineBuilder::FromPoint(Point pstart, const vector<double>& fun) const{
	//find value
	double value = approx->Val(pstart, fun);
	//assemble all edges of contour
	HMCont2D::PCollection points;
	HMCont2D::ECollection edges;
	for (int i=0; i<grid->n_cells(); ++i)
		AddLine(grid->get_cell(i), edges, points, value, fun);
	//collect into contours
	auto et = HMCont2D::Assembler::ETree(edges);
	//choose contour which passes pstart
	HMCont2D::Contour* cnt;
	if (et.cont_count() == 1) cnt = et.get_contour(0);
	else{
		for (int i=0; i<et.cont_count(); ++i){
			cnt = et.get_contour(i);
			if (ISZERO(std::get<4>(cnt->coord_at(pstart)))) break;
		}
	}
	
	return HMCont2D::Container<HMCont2D::Contour>::DeepCopy(*cnt);
}

HMCallback::FunctionWithCallback<HMFem::TAuxGrid3> HMFem::AuxGrid3;
double TAuxGrid3::step_estimate(const HMCont2D::ContourTree& tree, int nrec){
	double area = HMCont2D::ContourTree::Area(tree);
	double triarea = area/2/nrec;
	//from a area of equilateral triangle
	return sqrt(2.3094017677*triarea);
}

void TAuxGrid3::input(const HMCont2D::ContourTree& _tree, const vector<HMCont2D::Contour>& _constraints){
	clear();
	//deep copy edges input to internal structure
	HMCont2D::ContourTree::DeepCopy(_tree, tree);
	constraints.reserve(_constraints.size());
	for (auto& c: _constraints){
		constraints.emplace_back(new HMCont2D::Contour());
		HMCont2D::Contour::DeepCopy(c, *constraints.back());
	}
	//fill mandatory points
	for (int i=0; i<constraints.size(); ++i){
		for (int j=0; j<constraints.size(); ++j) if (i!=j){
			mandatory_intersections(*constraints[i], *constraints[j]);
		}
		for (auto& c: tree.nodes) mandatory_intersections(*constraints[i], *c);
	}
}

void TAuxGrid3::mandatory_intersections(HMCont2D::Contour& c1, HMCont2D::Contour& c2){
	auto cres = HMCont2D::Algos::CrossAll(c1, c2);
	for (auto& it: cres){
		Point& p = std::get<1>(it);
		auto gp1 = c1.GuaranteePoint(p, pcol);
		auto gp2 = c2.GuaranteePoint(p, pcol);
		mandatory_points.insert(std::get<1>(gp1));
		mandatory_points.insert(std::get<1>(gp2));
	}
}

void TAuxGrid3::adopt_contour(HMCont2D::Contour& cont, double h,
		vector<vector<Point>>& lost){
	auto repart = HMCont2D::Algos::Coarsening(cont, mandatory_points, pcol, h,
			ZEROANGLE, 0.2);
	vector<Point*> ret;
	for (int i=0; i<repart.size(); ++i){
		auto& v = repart[i];
		ret.push_back(v[0]);
		if (v.size() > 1){
			lost.push_back(vector<Point>());
			for (int i=0; i<v.size(); ++i){
				lost.back().push_back(*v[i]);
			}
			lost.back().push_back(*repart[(i+1) % repart.size()][0]);
		}
	}
	cont = HMCont2D::Constructor::ContourFromPoints(ret, cont.is_closed());
}

void TAuxGrid3::adopt_boundary(HMCont2D::ContourTree& tree, double h,
		vector<vector<Point>>& lost){
	for (auto n: tree.nodes) adopt_contour(*n, h, lost);
	tree.ReloadEdges();
}

Grid43 HMFem::TAuxGrid3::_run(const HMCont2D::ContourTree& _tree,
		const vector<HMCont2D::Contour>& _constraints,
		int nrec, int nmax){
	assert(nrec<nmax);
	//1)
	callback->step_after(5, "Estimate step");
	double hest = step_estimate(_tree, nrec);
	//2)
	callback->step_after(10, "Adopt boundary");
	input(_tree, _constraints);
	vector<vector<Point>> lost_points;
	adopt_boundary(tree, CORRECTION_FACTOR*hest, lost_points);
	for (auto c: constraints) adopt_contour(*c, CORRECTION_FACTOR*hest, lost_points);
	//3)
	callback->step_after(60, "Triangulate");
	auto ret = Grid43::Build3(tree, constraints, 10*CORRECTION_FACTOR*hest);
	if (ret->n_points() > nmax) { clear(); throw std::runtime_error("Failed to build auxiliary triangle grid");}
	if (lost_points.size() == 0) return std::move(*ret);
	//4)
	callback->step_after(15, "Snapping");
	ret->AddSegments(lost_points);
	//5)
	callback->step_after(10, "Boundary triangles");
	ret = Grid43::Build3(ret.get());
	//6)
	assert(ret->check());
	clear();
	return std::move(*ret);
}
void HMFem::TAuxGrid3::clear(){
	//explicit clear because there is the reusable entry of TAuxGrid class
	tree = HMCont2D::ContourTree();
	constraints.clear();
	pcol.clear();
	mandatory_points.clear();
}

Grid43 HMFem::TAuxGrid3::_run(const HMCont2D::ContourTree& tree, int nrec, int nmax){
	return _run(tree, {}, nrec, nmax);
}
Grid43 HMFem::TAuxGrid3::_run(const HMCont2D::Contour& cont, int nrec, int nmax){
	HMCont2D::ContourTree tree;
	tree.AddContour(cont);
	return _run(tree, nrec, nmax);
}
