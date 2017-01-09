#include "hmblay.hpp"
#include "procgrid.h"
#include "constructor.hpp"

GridGeom HalfCirc(HM2D::EdgeData& b1, HM2D::EdgeData& b2, double sz){
	Point b1f = *HM2D::Contour::First(b1);
	Point b1l = *HM2D::Contour::Last(b1);
	Point b2f = *HM2D::Contour::First(b2);
	Point b2l = *HM2D::Contour::Last(b2);
	//calculate size
	bool has_center = b1f == b2f;
	double rad = Point::dist(b1l, b2l)/2.0;
	int narc = 2*M_PI*rad/sz;
	if (narc<8) narc = 8;
	if (narc % 2 != 0) narc+=1;
	sz = 2*M_PI*rad/narc; 
	
	int nrad = b1.size();
	//if (!has_center) ++nrad;
	Point pc;
	if (has_center) pc = b1f;
	else pc = (b1f + b2f)/2.0;

	//build circular grid
	GridGeom ret = GGeom::Constructor::Circle(Point(0, 0), rad, narc, nrad, true);

	//change radii
	std::vector<double> rads;
	auto pts1 = HM2D::Contour::OrderedPoints(b1);
	auto pts2 = HM2D::Contour::OrderedPoints(b2);
	if (!has_center) pts1[0]->set(pc);
	if (!has_center) pts2[0]->set(pc);
	for (int i=1; i<pts1.size(); ++i){
		rads.push_back(Point::dist(pc, *pts1[i]));
	}
	for (int j=0; j<nrad; ++j){
		for (int i=0; i<narc; ++i){
			int k = (nrad-1-j)*narc + i;
			vecSetLen(*ret.get_point(k), rads[j]);
		}
	}

	//cut
	vector<const Cell*> rmcells;
	for (int j=0; j<nrad; ++j){
		for (int i=narc/2; i<narc; ++i){
			int k = (nrad-1-j)*narc+i;
			rmcells.push_back(ret.get_cell(k));
		}
	}
	GGeom::Modify::RemoveCells(ret, rmcells);
	
	//translate and rotate
	double angle = atan2(b1l.y - b1f.y, b1l.x - b1f.x);
	GGeom::Modify::PointModify(ret, [pc, angle](GridPoint* p){
				p->set(vecRotate(*p, angle));
				*p += pc;
			});
	//snap to input boundaries
	for (int j=0; j<nrad; ++j){
		int k1 = (nrad-1-j)*(narc/2+1) + 0;
		int k2 = (nrad-1-j)*(narc/2+1) + narc/2;
		ret.get_point(k1)->set(*pts1[j+1]);
		ret.get_point(k2)->set(*pts2[j+1]);
	}

	//connect to parent
	if (!has_center){
		//move
		Vect mv = b1l - b1f;
		mv = vecRotate(mv, M_PI/2.0);
		vecSetLen(mv, sz);
		GGeom::Modify::PointModify(ret, [mv](GridPoint* p){ *p += mv; });

		//additional section
		ShpVector<GridPoint> apoints, ap1, ap2;
		ShpVector<Cell> acells;
		for (auto p: pts1) aa::add_shared(ap1, GridPoint(*p));
		for (auto p: pts2) aa::add_shared(ap2, GridPoint(*p));
		ap1[0]->set(b1f); ap2[0]->set(b2f);
		std::copy(ap1.begin(), ap1.end(), std::back_inserter(apoints));
		std::copy(ap2.begin(), ap2.end(), std::back_inserter(apoints));
		for (int j=0; j<nrad; ++j){
			shared_ptr<Cell> c1(new Cell);
			int k1 = (nrad-j)*(narc/2+1) + 0;
			int k2 = (nrad-j-1)*(narc/2+1) + 0;
			c1->points.push_back(ap1[j].get());
			c1->points.push_back(ap1[j+1].get());
			c1->points.push_back(ret.get_point(k2));
			c1->points.push_back(ret.get_point(k1));
			acells.push_back(c1);
		}
		for (int j=0; j<nrad; ++j){
			shared_ptr<Cell> c1(new Cell);
			int k1 = (nrad-j)*(narc/2+1);
			int k2 = (nrad-j-1)*(narc/2+1) + narc/2;
			if (j!=0) k1 += narc/2;
			c1->points.push_back(ap2[j+1].get());
			c1->points.push_back(ap2[j].get());
			c1->points.push_back(ret.get_point(k1));
			c1->points.push_back(ret.get_point(k2));
			acells.push_back(c1);
		}
		//center triangle
		Cell* c = aa::add_shared(acells, Cell());
		c->points.push_back(ret.get_point(nrad*(narc/2+1)));
		c->points.push_back(ap2[0].get());
		c->points.push_back(ap1[0].get());
		GridGeom ag = GGeom::Constructor::FromData(apoints, acells);
		GGeom::Modify::ShallowAdd(&ag, &ret);

	}
	return ret;
}

HMCallback::FunctionWithCallback<HMBlay::TBuildStripeGrid> HMBlay::BuildStripeGrid;

GridGeom HMBlay::TBuildStripeGrid::_run(const HM2D::EdgeData& cont,
		const std::vector<double>& partition,
		int tip_algo, Point& bl, Point& br, Point& tr, Point& tl){
	HM2D::EdgeData cp(cont);
	HMBlay::Input opt;
	opt.partition = partition;
	if (opt.partition[0] > 0) opt.partition.insert(opt.partition.begin(), 0);
	opt.bnd_step_method = HMBlay::BndStepMethod::NO_BND_STEPPING;
	opt.edges = &cp;
	opt.start = *HM2D::Contour::First(cont);
	opt.end = *HM2D::Contour::Last(cont);

	callback->step_after(45, "Upper grid");
	opt.direction = HMBlay::Direction::INNER;
	GridGeom g1 = HMBlay::BuildBLayerGrid({opt});
	HM2D::EdgeData cleft1, cright1;
	{
		vector<Point> _pv1, _pv2;
		if (HM2D::Contour::IsOpen(cont)){
			int Nplen = cont.size()+1;
			for (int i=0; i<opt.partition.size(); ++i){
				_pv1.push_back(*g1.get_point(i*Nplen));
				_pv2.push_back(*g1.get_point(i*Nplen+Nplen-1));
			}
			cleft1 = HM2D::Contour::Constructor::FromPoints(_pv1);
			cright1 = HM2D::Contour::Constructor::FromPoints(_pv2);
			tl = *HM2D::Contour::Last(cleft1);
			tr = *HM2D::Contour::Last(cright1);
		} else { tl = *g1.get_point(g1.n_points()-1); tr = tl; }
	}

	callback->step_after(45, "Lower grid");
	opt.direction = HMBlay::Direction::OUTER;
	GridGeom g2 = HMBlay::BuildBLayerGrid({opt});
	HM2D::EdgeData cleft2, cright2;
	{
		vector<Point> _pv1, _pv2;
		if (HM2D::Contour::IsOpen(cont)){
			int Nplen = cont.size()+1;
			for (int i=0; i<opt.partition.size(); ++i){
				_pv2.push_back(*g2.get_point(i*Nplen));
				_pv1.push_back(*g2.get_point(i*Nplen+Nplen-1));
			}
			cleft2 = HM2D::Contour::Constructor::FromPoints(_pv1);
			cright2 = HM2D::Contour::Constructor::FromPoints(_pv2);
			bl = *HM2D::Contour::Last(cleft2);
			br = *HM2D::Contour::Last(cright2);
		} else { bl = *g2.get_point(g2.n_points()-1); br=bl; }
	}

	//if no center line
	if (partition[0] > 0){
		int N = cont.size(); if (HM2D::Contour::IsOpen(cont)) ++N;

		vector<int> lowpts1, lowpts2, lowcls2;
		for (int i=0; i<N; ++i) lowpts1.push_back(i);
		for (int i=0; i<N; ++i) lowpts2.push_back(i+N);
		if (HM2D::Contour::IsOpen(cont)) std::reverse(lowpts2.begin(), lowpts2.end());
		else std::reverse(lowpts2.begin()+1, lowpts2.end());
		for (int i=0; i<cont.size(); ++i) lowcls2.push_back(i);

		for (int i=0; i<N; ++i){
			g1.get_point(lowpts1[i])->set(*g2.get_point(lowpts2[i]));
		}

		vector<const Cell*> rc;
		for (auto i: lowcls2) rc.push_back(g2.get_cell(i));
		GGeom::Modify::RemoveCells(g2, rc);

		if (HM2D::Contour::IsOpen(cont)){
			cleft1.erase(cleft1.begin());
			cleft2.erase(cleft2.begin());
			cright1.erase(cright1.begin());
			cright2.erase(cright2.begin());
		}
		
	}

	callback->step_after(5, "Tip grids");
	shared_ptr<GridGeom> tip1, tip2;

	if (HM2D::Contour::IsClosed(cont)){ 
	} else if (tip_algo == 0){
	} else if (tip_algo == 1){
		if (cleft1.size()>0){
			tip1.reset(new GridGeom(HalfCirc(cleft1, cleft2, cont[0]->length())));
			tip2.reset(new GridGeom(HalfCirc(cright2, cright1, cont.back()->length())));
		} else {
			ShpVector<GridPoint> pts;
			ShpVector<Cell> cls;
			Vect av1 = tl - *HM2D::Contour::First(cont);
			vecSetLen(av1, (vecLen(av1) + cont[0]->length())/2.0);
			av1 = vecRotate(av1, M_PI/2.);
			aa::add_shared(pts, GridPoint(tl));
			aa::add_shared(pts, GridPoint( *HM2D::Contour::First(cont) + av1  ));
			aa::add_shared(pts, GridPoint(bl));
			aa::add_shared(cls, Cell());
			cls[0]->points = {pts[0].get(), pts[1].get(), pts[2].get()};
			tip1.reset(new GridGeom(GGeom::Constructor::FromData(pts, cls)));

			pts.clear(); cls.clear();
			Vect av2 = br - *HM2D::Contour::Last(cont);
			vecSetLen(av2, (vecLen(av2) + cont.back()->length())/2.0);
			av2 = vecRotate(av2, M_PI/2.);
			aa::add_shared(pts, GridPoint(br));
			aa::add_shared(pts, GridPoint( *HM2D::Contour::Last(cont) + av2  ));
			aa::add_shared(pts, GridPoint(tr));
			aa::add_shared(cls, Cell());
			cls[0]->points = {pts[0].get(), pts[1].get(), pts[2].get()};
			tip2.reset(new GridGeom(GGeom::Constructor::FromData(pts, cls)));
		}
	}

	callback->step_after(5, "Merge");
	GridGeom ret = GGeom::Constructor::EmptyGrid();
	GGeom::Modify::ShallowAdd(&g1, &ret);
	GGeom::Modify::ShallowAdd(&g2, &ret);
	if (tip1) GGeom::Modify::ShallowAdd(tip1.get(), &ret);
	if (tip2) GGeom::Modify::ShallowAdd(tip2.get(), &ret);
	GGeom::Repair::Heal(ret);

	return ret;
}
