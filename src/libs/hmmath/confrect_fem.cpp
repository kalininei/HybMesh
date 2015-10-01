#include "confrect_fem.hpp"
#include "hmfem.hpp"
#include "femassembly.hpp"

using namespace HMMath::Conformal::Impl::ConfFem;

ToRect::ToRect(const vector<Point>& path, int i1, int i2, int i3, double h){
	//1) build fem grid: fill grid, origs, HMCont2D::Contour data
	BuildGrid(path, i1, i2, i3, h);
	//2) compute mapping: fill u, v, module
	DoMapping();
	//3) build inverse
	BuildInverse();
}

#include "fileproc.h"
void ToRect::BuildGrid(const vector<Point>& path, int i1, int i2, int i3, double h){
	grid = HMFem::Grid43::Build3(path, h);
	approx = grid->GetApprox();
	save_vtk(*grid, "_dbgout.vtk");
	//2) restore origs
	HMCont2D::ContourTree tree = GGeom::Info::Contour(*grid);
	assert(tree.cont_count() == 1);
	auto& c = tree.nodes[0];
	auto bp = c->ordered_points();
	for (auto& p: path){
		auto fnd = std::find_if(bp.begin(), bp.end(),
			[&p](const Point* sp){ return *sp == p; });
		assert(fnd != bp.end());
		origs.push_back(static_cast<const GridPoint*>(*fnd)->get_ind());
	}
	//3) build orig contours
	left = HMCont2D::Contour::Assemble(*c, grid->get_point(origs[0]), grid->get_point(origs[i1]));
	bottom = HMCont2D::Contour::Assemble(*c, grid->get_point(origs[i1]), grid->get_point(origs[i2]));
	right = HMCont2D::Contour::Assemble(*c, grid->get_point(origs[i2]), grid->get_point(origs[i3]));
	top = HMCont2D::Contour::Assemble(*c, grid->get_point(origs[i3]), grid->get_point(origs[0]));
}

std::array<vector<const GridPoint*>, 4> ToRect::BndPoints() const{
	auto bpvec = GGeom::Info::BoundaryPoints(*grid);
	std::list<const GridPoint*> bplist;
	for (auto& x: bpvec) bplist.push_back(x.get());

	//remove original points from list
	std::remove_if(bplist.begin(), bplist.end(), [&](const GridPoint* p){
		return std::find(origs.begin(), origs.end(), p->get_ind())
			!= origs.end();
	});

	//assemble result
	std::array<vector<const GridPoint*>, 4> ret;
	auto& rleft=ret[0], &rbot=ret[1], &rright=ret[2], &rtop=ret[3];

	//fills v from points in bplist and orig points according to their position in cont.
	//points are removed from bplist after their addition to v.
	auto gather = [&](const HMCont2D::Contour& cont, vector<const GridPoint*>& v){
		//1) sort out points from blist according to bounding box
		BoundingBox bbox = HMCont2D::Contour::BBox(cont);
		std::list<const GridPoint*> so = bbox.Filter(bplist, 1, 1, 0);
		//2) for each contour edge find points
		auto bp = cont.ordered_points();
		std::map<double, const GridPoint*> gp;
		for (int i=0; i<bp.size()-1; ++i){
			Vect v1= *bp[i+1] - *bp[i]; vecNormalize(v1);
			double Dor=v1.y*bp[i]->x-v1.x*bp[i]->y;
			bool byx = ISEQ(bp[i]->x, bp[i+1]->x);
			double cl = (byx) ? (bp[i+1]->y - bp[i]->y) : bp[i+1]->x - bp[i]->x;
			std::remove_if(so.begin(), so.end(), [&](const GridPoint* p){
				//define ksi: wheit of point on section
				double ksi;
				double Dp =v1.y*p->x-v1.x*p->y;
				if (fabs(Dor-Dp)>geps) ksi = 10;
				else{
					if (byx) ksi = (p->y-bp[i]->y)/cl;
					else     ksi = (p->x-bp[i]->x)/cl;
				}
				//process
				if (ksi>0 && ksi<1){
					bplist.remove(p);
					gp[(double)i+ksi] = p;
					return true;
				} else return false;
			});
		}
		//3) assemble
		auto it = gp.begin();
		for (int i=0; i<bp.size()-1; ++i){
			v.push_back(static_cast<const GridPoint*>(bp[i]));
			while (it!=gp.end() && it->first < i+1)
				v.push_back((it++)->second);
		}
		v.push_back(static_cast<const GridPoint*>(bp.back()));
	};

	gather(left, rleft);
	gather(right, rright);
	gather(bottom, rbot);
	gather(top, rtop);

	return ret;
}

void ToRect::DoMapping(){
	using namespace HMFem;
	auto bp = BndPoints();
	u.resize(grid->n_points(), 0.0);
	v.resize(grid->n_points(), 0.0);
	//1) assemble pure laplas operator
	auto laplas = Assemble::PureLaplas(*grid);
	//2.1) Laplas problem for u
	auto ulaplas = LaplasProblem(grid, laplas);
	ulaplas.SetDirichlet(bp[0], [](const GridPoint*){ return 0; });
	ulaplas.SetDirichlet(bp[2], [](const GridPoint*){ return 1; });
	ulaplas.Solve(u);
	//2.2) Laplas problem for v
	auto vlaplas = LaplasProblem(grid, laplas);
	vlaplas.SetDirichlet(bp[1], [](const GridPoint*){ return 0; });
	vlaplas.SetDirichlet(bp[3], [](const GridPoint*){ return 1; });
	vlaplas.Solve(v);

	//3) Compute module
	_module = -vlaplas.IntegralDfDn(bp[1], v);

	//4) correct u according to module
	for (auto& x: u) x *= _module;
}

void ToRect::BuildInverse(){
	//grid
	inv_grid.reset(new HMFem::Grid43());
	GGeom::Modify::DeepAdd(grid.get(), inv_grid.get());
	auto modfun = [&](GridPoint* p){
		p->x = u[p->get_ind()];
		p->y = v[p->get_ind()];
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

double ToRect::HEstimate(const vector<Point>& path, int segn, int nmax){
	double h = 1e10;
	//1) calculate h depending on segment partition
	for (int i=0; i<path.size(); ++i){
		int inext = (i == path.size() - 1) ? 0 : i+1;
		double d = Point::dist(path[i], path[inext]);
		if (d<h) h = d;
	}
	//using +1 because of TriGrid::TriangulateArea implementation
	h /= (segn+1);
	//2) estimate total number of cells and correct h if necessary
	double A = 0;
	for (int i=1; i<path.size()-1; ++i){
		A += triarea(path[0], path[i], path[i+1]);
	}
	A = fabs(A);
	int Nest = A/(h*h);
	if (Nest<=nmax) return h;
	else return sqrt(A/nmax);
}

shared_ptr<ToRect>
ToRect::Build(const vector<Point>& path, int i1, int i2, int i3, const Options& opt){
	// - calculate step size
	double h = HEstimate(path, opt.fem_segment_partition, opt.fem_nmax);
	// - build mapping
	shared_ptr<ToRect> ret(new ToRect(path, i1, i2, i3, h));
	// - return
	if (ret->module() > 0) return ret;
	else return 0;
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





