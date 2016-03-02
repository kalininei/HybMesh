#include "femgrid43.hpp"
#include "trigrid.h"

using namespace HMFem;

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

	Point ksieta;
	ksieta.x = ( j22*(p.x - x1) - j21*(p.y - y1))/modj;
	ksieta.y = (-j12*(p.x - x1) + j11*(p.y - y1))/modj;

	assert( ksieta.x > -0.1 && ksieta.x < 1.1 );
	assert( ksieta.y > -0.1 && ksieta.y < 1-ksieta.x+0.1 );

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

double Grid43::Approximator::Val(Point p, const vector<double>& fun) const{
	const Cell* c = cfinder->Find(p);
	Point ksieta = LocalCoordinates(c, p);
	return Interpolate(c, ksieta, fun);
}


vector<double> Grid43::Approximator::Vals(Point p, const vector<const vector<double>*>& funs) const{
	const Cell* c = cfinder->Find(p);
	Point ksieta = LocalCoordinates(c, p);
	vector<double> ret;
	for (auto& fun: funs) ret.push_back(Interpolate(c, ksieta, *fun));
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

