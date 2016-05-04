#include "confrect_fem.hpp"
#include "hmfem.hpp"
#include "femassembly.hpp"
#include "debug_grid2d.h"

#define FEM_NMAX 1000000
#define APPROX_PART 200

using namespace HMMath::Conformal::Impl::ConfFem;
namespace{

double IntegralGrad2(shared_ptr<HMFem::Grid43> grid, const vector<double>& v){
	using namespace HMFem;
	shared_ptr<Mat> ddx = Assemble::DDx(*grid);
	shared_ptr<Mat> ddy = Assemble::DDy(*grid);
	vector<double> vdx(ddx->rows()), vdy(ddx->rows()), mass = Assemble::LumpMass(*grid);
	ddx->MultVec(v, vdx); 
	ddy->MultVec(v, vdy);
	for (int i=0; i<mass.size(); ++i){ vdx[i] /= mass[i]; vdy[i] /= mass[i]; }
	//3.2) |grad|^2 -> vdx
	std::transform(vdx.begin(), vdx.end(), vdy.begin(), vdx.begin(),
			[](double a, double b){ return a*a + b*b; });
	//3.3) integral
	return std::inner_product(vdx.begin(), vdx.end(), mass.begin(), 0.0);
}

void GradVectors(shared_ptr<HMFem::Grid43> grid, const vector<double>& v,
		vector<double>& dx, vector<double>& dy){
	using namespace HMFem;
	shared_ptr<Mat> ddx = Assemble::DDx(*grid);
	shared_ptr<Mat> ddy = Assemble::DDy(*grid);
	vector<double> mass = Assemble::LumpMass(*grid);
	dx.resize(grid->n_points()); dy.resize(grid->n_points());
	ddx->MultVec(v, dx); 
	ddy->MultVec(v, dy);
	for (int i=0; i<mass.size(); ++i){ dx[i] /= mass[i]; dy[i] /= mass[i]; }
}

vector<int> GetGridIndicies(const vector<Point>& path, const vector<Point*>& cont){
	vector<int> origs;
	for (auto& p: path){
		auto fnd = std::find_if(cont.begin(), cont.end(),
			[&p](const Point* sp){return *sp == p;});
		assert(fnd != cont.end());
		origs.push_back(static_cast<const GridPoint*>(*fnd)->get_ind());
	}
	return origs;
}

vector<int> GetGridIndicies(const vector<Point*>& cont){
	vector<int> origs;
	for (auto& p: cont){
		origs.push_back(static_cast<const GridPoint*>(p)->get_ind());
	}
	return origs;
}

vector<const GridPoint*> ExtractPoints(const GridGeom& grid, const vector<int>& ind){
	vector<const GridPoint*> ret;
	for (int i: ind) ret.push_back(grid.get_point(i));
	return ret;
}

double LinTriangleSize(double area, int ncells){ return 1.5*sqrt(4.0*(area/ncells)/sqrt(3));}


}

HMCallback::FunctionWithCallback<ToRect::TBuild> ToRect::Build;
ToRect ToRect::TBuild::_run(const vector<Point>& path, int i1, int i2, int i3, const Options& opt){
	ToRect ret;

	//1) build fem grid: fill grid, origs, HMCont2D::Contour data
	BuildGrid(ret, path, i1, i2, i3, opt.fem_nrec);
	//2) compute mapping: fill u, v, module, approx
	DoMapping(ret);
	//3) build inverse
	BuildInverse(ret);

	return ret;
}

void ToRect::TBuild::BuildGrid(ToRect& r, const vector<Point>& path, int i1, int i2, int i3, int n){
	//1) build grid
	auto subcaller = callback->bottom_line_subrange(45);
	auto cont = HMCont2D::Constructor::ContourFromPoints(path, true);
	r.grid.reset(new HMFem::Grid43(HMFem::AuxGrid3.UseCallback(subcaller, cont, n, FEM_NMAX)));

	callback->step_after(5, "Assembling original contour");
	r.approx = r.grid->GetApprox();
	//2) restore origs
	HMCont2D::Contour c = GGeom::Info::Contour1(*r.grid);
	r.origs = GetGridIndicies(path, c.ordered_points());
	//3) build orig contours
	auto left = HMCont2D::Assembler::Contour1(c, r.grid->get_point(r.origs[0]), r.grid->get_point(r.origs[i1]));
	auto bottom = HMCont2D::Assembler::Contour1(c, r.grid->get_point(r.origs[i1]), r.grid->get_point(r.origs[i2]));
	auto right = HMCont2D::Assembler::Contour1(c, r.grid->get_point(r.origs[i2]), r.grid->get_point(r.origs[i3]));
	auto top = HMCont2D::Assembler::Contour1(c, r.grid->get_point(r.origs[i3]), r.grid->get_point(r.origs[0]));
	r.ileft = GetGridIndicies(left.ordered_points());
	r.iright = GetGridIndicies(right.ordered_points());
	r.itop = GetGridIndicies(top.ordered_points());
	r.ibottom = GetGridIndicies(bottom.ordered_points());
}

void ToRect::TBuild::DoMapping(ToRect& r){
	using namespace HMFem;
	r.u.resize(r.grid->n_points(), 0.0);
	r.v.resize(r.grid->n_points(), 0.0);
	//1) assemble pure laplas operator
	callback->step_after(10, "Laplace operator");
	auto laplas = Assemble::PureLaplas(*r.grid);
	//2.1) Laplas problem for u
	callback->step_after(10, "U-problem");
	auto ulaplas = LaplasProblem(r.grid, laplas);
	ulaplas.SetDirichlet(ExtractPoints(*r.grid, r.ileft), [](const GridPoint*){ return 0; });
	ulaplas.SetDirichlet(ExtractPoints(*r.grid, r.iright), [](const GridPoint*){ return 1; });
	ulaplas.Solve(r.u);
	//2.2) Laplas problem for v
	callback->step_after(10, "V-problem");
	auto vlaplas = LaplasProblem(r.grid, laplas);
	vlaplas.SetDirichlet(ExtractPoints(*r.grid, r.ibottom), [](const GridPoint*){ return 0; });
	vlaplas.SetDirichlet(ExtractPoints(*r.grid, r.itop), [](const GridPoint*){ return 1; });
	vlaplas.Solve(r.v);

	//3) Compute modulus
	callback->step_after(5, "Modulus");
	//Calculating using area integration.
	//Taking into account 
	//    Int(dvdn)dtop = Int (dvdn*v) dG = Int (grad v)^2 dD
	r._module = IntegralGrad2(r.grid, r.v);
	
	//4) correct u according to modulus
	for (auto& x: r.u) x *= r._module;
}

void ToRect::TBuild::BuildInverse(ToRect& r){
	callback->step_after(15, "Inversion");
	//grid
	r.inv_grid.reset(new HMFem::Grid43());
	GGeom::Modify::DeepAdd(r.grid.get(), r.inv_grid.get());
	auto modfun = [&r](GridPoint* p){
		p->x = r.u[p->get_ind()];
		p->y = r.v[p->get_ind()];
	};
	GGeom::Modify::PointModify(*r.inv_grid, modfun);
	//approximator
	r.inv_approx = r.inv_grid->GetApprox();
	//inverse functions
	r.inv_u.resize(r.inv_grid->n_points());
	r.inv_v.resize(r.inv_grid->n_points());
	for (int i=0; i<r.inv_u.size(); ++i){
		r.inv_u[i] = r.grid->get_point(i)->x;
		r.inv_v[i] = r.grid->get_point(i)->y;
	}
}

vector<Point> ToRect::MapToPolygon(const vector<Point>& input) const{
	vector<Point> ret;
	vector<const vector<double>*> funs {&inv_u, &inv_v};
	for (auto& p: input){
		auto s = inv_approx->Vals(p, funs);
		ret.push_back(Point(s[0], s[1]));
	}
	return ret;
}

vector<Point> ToRect::MapToRectangle(const vector<Point>& input) const{
	vector<Point> ret;
	vector<const vector<double>*> funs {&u, &v};
	for (auto& p: input){
		auto s = approx->Vals(p, funs);
		ret.push_back(Point(s[0], s[1]));
	}
	return ret;
}

vector<Point> ToRect::RectPoints() const{
	vector<Point> ret;
	for (auto i: origs) ret.push_back(Point(u[i], v[i]));
	return ret;
}


// ================== Annulus
shared_ptr<ToAnnulus>
ToAnnulus::Build(const vector<Point>& outer_path, const vector<Point>& inner_path, const Options& opt){
	// - build mapping
	shared_ptr<ToAnnulus> ret(new ToAnnulus(outer_path, inner_path, opt));
	// - return
	if (ret->module() > 0) return ret;
	else return 0;
}

void ToAnnulus::BuildGrid1(const vector<Point>& outer_path, const vector<Point>& inner_path,
		int n){
	//grid
	auto c1 = HMCont2D::Constructor::ContourFromPoints(inner_path, true);
	auto c2 = HMCont2D::Constructor::ContourFromPoints(outer_path, true);
	HMCont2D::ContourTree pretree;
	pretree.AddContour(c1);
	pretree.AddContour(c2);
	grid.reset(new HMFem::Grid43(HMFem::AuxGrid3(pretree, n, FEM_NMAX)));
	approx = grid->GetApprox();
	//boundary indicies
	HMCont2D::ContourTree ct = GGeom::Info::Contour(*grid);
	auto outerc = ct.roots()[0]->ordered_points();
	auto innerc  = ct.roots()[0]->children[0]->ordered_points();
	//all bnd points
	outer = GetGridIndicies(outerc);
	inner = GetGridIndicies(innerc);
}

void ToAnnulus::DoMappingU(){
	//laplas problem
	auto ulaplas = HMFem::LaplasProblem(grid);
	//dirichlet: one on inner, zero on outer
	ulaplas.SetDirichlet(ExtractPoints(*grid, outer), [](const GridPoint*){return 0.0;});
	ulaplas.SetDirichlet(ExtractPoints(*grid, inner), [](const GridPoint*){return 1.0;});
	//solution
	u.resize(grid->n_points(), 0.0);
	ulaplas.Solve(u);
	//Modulus
	_module = exp( -2*M_PI/IntegralGrad2(grid, u) );
}

HMCont2D::Container<HMCont2D::Contour> ToAnnulus::SteepestDescentUCurve(){
	//get derivatives
	vector<double> dx, dy;
	GradVectors(grid, u, dx, dy);
	//step
	double h = sqrt(grid->area()/grid->n_cells())/2.0;
	//outer contour
	HMCont2D::Contour outerc = *GGeom::Info::Contour(*grid).roots()[0];
	//go
	vector<Point> ret; ret.push_back(pzero);
	while (outerc.IsWithin(ret.back())){
		vector<double> grad = approx->Vals(ret.back(), {&dx, &dy});
		Vect gradv {grad[0], grad[1]};
		vecSetLen(gradv, h);
		ret.push_back(ret.back() - gradv);
	}
	//assemble contour
	auto cont = HMCont2D::Constructor::ContourFromPoints(ret);
	//if last point lies outside outerc cut cont
	if (outerc.IsWithout(ret.back())){
		auto cr = HMCont2D::Algos::Cross(outerc, cont);
		*(cont.last()) = std::get<1>(cr);
	}
	//if distance beween last point and one before last is very short
	//remove point before last
	if (cont.data.back()->length()<h/4){
		Point *pp = cont.last();
		cont.data.pop_back();
		*cont.last() = *pp;
	}
	return cont;
}


void ToAnnulus::BuildGrid2(const vector<Point>& outer_path,
		const vector<Point>& inner_path,
		int n,
		const HMCont2D::Contour& razor){
	//starting from grid/approx filled with doubly-connected grid
	//0) assemble inner and outer contours
	auto innerc = std::make_shared<HMCont2D::Container<HMCont2D::Contour>>(
			HMCont2D::Constructor::ContourFromPoints(inner_path, true));
	auto outerc = std::make_shared<HMCont2D::Container<HMCont2D::Contour>>(
			HMCont2D::Constructor::ContourFromPoints(outer_path, true));
	HMCont2D::ContourTree ct;
	shared_ptr<HMCont2D::Contour> shi(innerc);
	shared_ptr<HMCont2D::Contour> sho(outerc);
	ct.AddContour(shi);
	ct.AddContour(sho);
	
	//1) assemble doubly-connected contour
	grid.reset(new HMFem::Grid43(HMFem::AuxGrid3(
		ct, vector<HMCont2D::Contour> {razor}, n, FEM_NMAX)));
	//2) rip double-connected contour by razor polyline
	//2.1) get points which lie on razor
	std::vector<const GridPoint*> razor_points;
	auto modfun1 = [&](GridPoint* p){
		auto ca = razor.coord_at(*p);
		if (ISZERO(std::get<4>(ca))) razor_points.push_back(p);
	};
	GGeom::Modify::PointModify(*grid, modfun1);
	//2.2) double razor points and add em to grid
	std::map<const GridPoint*, GridPoint*> m12;
	ShpVector<GridPoint> pcol;
	for (auto p: razor_points){
		shared_ptr<GridPoint> np(new GridPoint(*p));
		pcol.push_back(np);
		m12[p] = np.get();
	}
	GGeom::Modify::ShallowAdd(pcol, grid.get());

	//2.3) change razor->razor2 for cells which lie to the right
	auto modfun2 = [&](Cell* c){
		auto fnd0 = m12.find(c->get_point(0));
		auto fnd1 = m12.find(c->get_point(1));
		auto fnd2 = m12.find(c->get_point(2));
		if (fnd0 == m12.end() && fnd1 == m12.end() && fnd2 == m12.end()) return;
		Point cnt = (*c->get_point(0) + *c->get_point(1) + *c->get_point(2)) / 3.0;
		auto ca = razor.coord_at(cnt);
		HMCont2D::Edge red = *razor.edge(std::get<2>(ca));
		if (!razor.correctly_directed_edge(std::get<2>(ca))) red.Reverse();
		double ar = triarea(*red.pstart, *red.pend, cnt);
		if (ar > 0) return; //if cell is on the right side don't change
		if (fnd0 != m12.end()) c->points[0] = fnd0->second;
		if (fnd1 != m12.end()) c->points[1] = fnd1->second;
		if (fnd2 != m12.end()) c->points[2] = fnd2->second;
	};
	GGeom::Modify::CellModify(*grid, modfun2);
	//3) origs
	HMCont2D::Contour ct2 = GGeom::Info::Contour1(*grid);
	vector<Point*> allp = ct2.ordered_points(); allp.pop_back();
	outer_origs = GetGridIndicies(outer_path, allp);
	inner_origs = GetGridIndicies(inner_path, allp);
	//4) inner/outer
	inner.clear(); outer.clear();
	for (auto p: allp){
		if (ISZERO(std::get<4>(outerc->coord_at(*p)))){
			outer.push_back(static_cast<GridPoint*>(p)->get_ind());
		}
		if (ISZERO(std::get<4>(innerc->coord_at(*p)))){
			inner.push_back(static_cast<GridPoint*>(p)->get_ind());
		}
	}
	//5) razor io/oi
	raz_io.clear(); raz_oi.clear();
	for (auto m: m12){
		raz_io.push_back(m.first->get_ind());
		raz_oi.push_back(m.second->get_ind());
	}
	//6) approx
	approx = grid->GetApprox();
}

void ToAnnulus::DoMapping(){
	//laplas matrix
	auto lap = HMFem::Assemble::PureLaplas(*grid);
	//laplas problem for v
	auto vlaplas = HMFem::LaplasProblem(grid, lap);
	//dirichlet: 2*pi on razor_oi, zero on razor_io
	vlaplas.SetDirichlet(ExtractPoints(*grid, raz_io), [](const GridPoint*){return 0.0;});
	vlaplas.SetDirichlet(ExtractPoints(*grid, raz_oi), [](const GridPoint*){return 2.0*M_PI;});
	//solution
	v.resize(grid->n_points(), 0.0);
	vlaplas.Solve(v);

	//laplas problem for u
	auto ulaplas = HMFem::LaplasProblem(grid, lap);
	//dirichlet: 1 on top, _module on bot
	ulaplas.SetDirichlet(ExtractPoints(*grid, outer), [](const GridPoint*){return 1;});
	ulaplas.SetDirichlet(ExtractPoints(*grid, inner), [&](const GridPoint*){return _module;});
	//solution
	u.resize(grid->n_points(), 0.0);
	ulaplas.Solve(u);
}

void ToAnnulus::BuildInverse(){
	//ugrid
	inv_grid.reset(new HMFem::Grid43());
	GGeom::Modify::DeepAdd(grid.get(), inv_grid.get());
	auto modfun = [&](GridPoint* p){
		double &r = u[p->get_ind()], &phi = v[p->get_ind()];
		p->x = r*cos(phi);
		p->y = r*sin(phi);
	};
	GGeom::Modify::PointModify(*inv_grid, modfun);
	//approximator
	inv_approx = inv_grid->GetApprox();
	//inverse functions
	inv_u.resize(inv_grid->n_points());
	inv_v.resize(inv_grid->n_points());
	for (int i=0; i<inv_u.size(); ++i){
		inv_u[i] = grid->get_point(i)->x;
		inv_v[i] = grid->get_point(i)->y;
	}
}

ToAnnulus::ToAnnulus(const vector<Point>& outer_path, const vector<Point>& inner_path,
		const Options& opt){
	pzero = inner_path[0];
	//1) Build a grid for modulus+cirlce computation
	BuildGrid1(outer_path, inner_path, opt.fem_nrec);
	//2) Compute u and modulus
	DoMappingU();
	//3) Build a steepest descent curve
	auto sd = SteepestDescentUCurve();
	//4) Build a grid for conjugate function solution
	BuildGrid2(outer_path, inner_path, opt.fem_nrec, sd);
	//5) Recompute u on singly connected grid, Compute v
	DoMapping();
	//7) Build inverse grid
	BuildInverse();
}

vector<Point> ToAnnulus::InnerCirclePoints() const{
	vector<Point> ret;
	for (int i=0; i<inner_origs.size(); ++i){
		int j = inner_origs[i];
		ret.push_back(Point(u[j], v[j]));
	}
	return ret;
}

vector<Point> ToAnnulus::OuterCirclePoints() const{
	vector<Point> ret;
	for (int i=0; i<outer_origs.size(); ++i){
		int j = outer_origs[i];
		ret.push_back(Point(u[j], v[j]));
	}
	return ret;
}

vector<Point> ToAnnulus::MapToAnnulus(const vector<Point>& input) const{
	vector<Point> ret;
	vector<const vector<double>*> funs {&u, &v};
	for (auto& p: input){
		auto s = approx->Vals(p, funs);
		double x = s[0]*cos(s[1]);
		double y = s[0]*sin(s[1]);
		ret.push_back(Point(x, y));
	}
	return ret;
}

vector<Point> ToAnnulus::MapToOriginal(const vector<Point>& input) const{
	vector<Point> ret;
	vector<const vector<double>*> funs {&inv_u, &inv_v};
	vector<double> s(2);
	for (auto& p: input){
		try{
			s = inv_approx->Vals(p, funs);
		} catch (GGeom::EOutOfArea& e){
			//if point was not found due to circle geom. approximation error
			//project point to inv_grid boundary and try again
			double rad = vecLen(p);
			Point pnew;
			if (fabs(rad-1.0)<geps || fabs(rad-_module)<geps){
				pnew = HMCont2D::ECollection::ClosestPoint(*InvGridContour(), p);
				s = inv_approx->Vals(pnew, funs);
			} else throw;
		}
		ret.push_back(Point(s[0], s[1]));
	}
	return ret;
}

double ToAnnulus::PhiInner(int i) const{
	return v[inner_origs[i]];
}
double ToAnnulus::PhiOuter(int i) const{
	return v[outer_origs[i]];
}

const HMCont2D::Contour* ToAnnulus::InvGridContour() const{
	if (!_inv_cont){
		_inv_cont = std::make_shared<HMCont2D::Contour>(
				GGeom::Info::Contour1(*inv_grid));
	}
	return _inv_cont.get();
}
