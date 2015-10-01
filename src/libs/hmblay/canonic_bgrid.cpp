#include "canonic_bgrid.hpp"
#include "hmfem.hpp"

using namespace HMBlay::Impl;

shared_ptr<MappedRect>
MappedRect::Factory(HMCont2D::Contour& left, HMCont2D::Contour& right,
		HMCont2D::Contour& bottom, HMCont2D::Contour& top,
		bool use_rect_approx){
	assert(*left.first() == *bottom.first());
	assert(*left.last() == *top.first());
	assert(*right.first() == *bottom.last());
	assert(*right.last() == *top.last());

	//if left/right coincide we need mapping to a ring
	auto leftcrn = left.corner_points();
	auto rightcrn = right.corner_points();
	if (leftcrn.size() == rightcrn.size()){
		if (std::equal(leftcrn.begin(), leftcrn.end(), rightcrn.begin(),
				[](Point* p1, Point* p2){ return *p1 == *p2; })){
			return RectForClosedArea::Build(bottom, top);
		}
	}
	//left/right do not coincide. We need mapping to rectangle.
	return std::make_shared<RectForOpenArea>(
		left, right, bottom, top, use_rect_approx
	);
};

shared_ptr<MappedRect>
MappedRect::Factory(HMCont2D::Contour& left, HMCont2D::Contour& right,
		HMCont2D::Contour& bottom, double h,
		bool use_rect_approx){
	assert(*left.first() == *bottom.first());
	assert(*right.first() == *bottom.last());

	//1) if bottom is a closed line draw outer contour at h from bottom
	if (*bottom.first() == *bottom.last()){
		return RectForClosedArea::Build(bottom, bottom.first(), h);
	}

	//2) cut left and right to size h since we need it
	auto left2 = HMCont2D::Contour::CutByLen(left, 0, h);
	auto right2 = HMCont2D::Contour::CutByLen(right, 0, h);
	bool straight_left = left2.is_straight();
	bool straight_right = right2.is_straight();
	bool straight_bottom = bottom.is_straight();

	//3) if all lines are straight and are at angles (pi/2 +/- pi/16)
	//   -> simply connect end points
	vector<Point*> op_left = left2.ordered_points();
	vector<Point*> op_right = right2.ordered_points();
	vector<Point*> op_bottom = bottom.ordered_points();
	double an1 = Angle(*op_left[1], *op_left[0], *op_bottom[1]);
	double an2 = Angle(*op_bottom[op_bottom.size()-2], *op_right[0], *op_right[1]);
	auto is90 = [](double a)->bool { double u = a*32/M_PI; return u>14.9 && u<17.1; };
	if (straight_left && straight_right && straight_bottom &&
			is90(an1), is90(an2)){
		HMCont2D::Contour top;
		top.add_value(HMCont2D::Edge(op_left.back(), op_right.back()));
		return Factory(left2, right2, bottom, top, use_rect_approx);
	}

	//4) if left and right are straight and angles are >= pi/2
	//   -> draw equidistant polygon around bot
	HMCont2D::Container<HMCont2D::Contour> top;
	if (straight_left && straight_right &&
			ISEQGREATER(an1, M_PI/2) && ISEQGREATER(an2, M_PI/2)){
		auto ellipse = HMCont2D::Contour::OffsetOuter(bottom, h);
		auto pl = ellipse.GuaranteePoint(*left2.last(), ellipse.pdata);
		auto pr = ellipse.GuaranteePoint(*right2.last(), ellipse.pdata);
		top = HMCont2D::Constructor::CutContour(ellipse, *right2.last(), *left2.last());
		top.ReallyReverse();
		*left2.last() = *top.first();
		*right2.last() = *top.last();
		goto GEOMETRY_RESULT_CHECK;
	}

	//6) try weighted geometric procedures
	_THROW_NOT_IMP_;

GEOMETRY_RESULT_CHECK:
	if (top.size() > 0){
		//no crosses for bottom,
		//cross at the end point for left, right
		auto cr1 = HMCont2D::Contour::Cross(bottom, top);
		if (std::get<0>(cr1) == true) goto GEOMETRY_FAILED;
		auto cr2 = HMCont2D::Contour::Cross(left2, top);
		if (std::get<0>(cr2) && !ISEQ(std::get<2>(cr2), 1.0)) goto GEOMETRY_FAILED;
		auto cr3 = HMCont2D::Contour::Cross(right2, top);
		if (std::get<0>(cr3) && !ISEQ(std::get<2>(cr3), 1.0)) goto GEOMETRY_FAILED;
		//all checks were passed
		return Factory(left2, right2, bottom, top, use_rect_approx);
	}

GEOMETRY_FAILED:
	//7) Failed to connect nodes geometrically.
	//   We need a really smart algorithm here. Smth like
	//   D. Gaier, On an area problem in conformal mapping, Results in Mathematics 10 (1986) 66-81.
	_THROW_NOT_IMP_;
}

double MappedRect::top2bot(double w) const{ return conf2bot(top2conf(w)); }
double MappedRect::bot2top(double w) const{ return conf2top(bot2conf(w)); }


// ============ Open area
RectForOpenArea::RectForOpenArea(HMCont2D::Contour& left, HMCont2D::Contour& right,
			HMCont2D::Contour& bottom, HMCont2D::Contour& top, bool use_rect_approx):
			MappedRect(left, right, bottom, top, use_rect_approx){
	//1. assemble polygon
	std::vector<Point> path;
	std::array<int, 4> corners;
	corners[0] = 0;
	auto tmp = left.corner_points(); std::reverse(tmp.begin(), tmp.end());
	for (Point* p: tmp) path.push_back(*p);
	path.pop_back();
	corners[1] = path.size();
	for (Point* p: bottom.corner_points()) path.push_back(*p);
	path.pop_back();
	corners[2] = path.size();
	for (Point* p: right.corner_points()) path.push_back(*p);
	path.pop_back();
	corners[3] = path.size();
	tmp = top.corner_points(); std::reverse(tmp.begin(), tmp.end());
	for (Point* p: tmp) path.push_back(*p);
	path.pop_back();
	//2. build conformal transformation into rectangle [0,m]x[0,1]
	HMMath::Conformal::Options opt;
	opt.use_rect_approx = use_rect_approx;
	opt.right_angle_eps = M_PI/8.0;
	opt.length_weight = 1.05;
	core = HMMath::Conformal::Rect::Factory(path, corners, opt);
}

HMCont2D::PCollection RectForOpenArea::MapToReal(const vector<const Point*>& p) const{
	vector<Point> ret(p.size());
	double m = core->module();
	//stretch to conformal rectangle
	std::transform(p.begin(), p.end(), ret.begin(), [&m](const Point* it){
			return Point(it->x*m, it->y);
	});
	//run transform
	ret = core->MapToPolygon(ret);
	//return
	HMCont2D::PCollection col;
	for (auto& it: ret){
		col.add_value(it);
	}
	return col;
}

HMCont2D::PCollection RectForOpenArea::MapToSquare(const vector<const Point*>& p) const{
	vector<Point> ret(p.size());
	std::transform(p.begin(), p.end(), ret.begin(), [](const Point* it){
			return Point(it->x, it->y);
	});
	//run transform
	ret = core->MapToRectangle(ret);
	//stretch to conformal rectangle
	double m = core->module();
	for (auto& p: ret) p.x/=m;
	//return
	HMCont2D::PCollection col;
	for (auto& it: ret){
		col.add_value(it);
	}
	return col;
}

double RectForOpenArea::conf2top(double w) const {
	auto p = MappedRect::MapToReal(Point(w, 1));
	return std::get<1>(TopContour().coord_at(p));
};
double RectForOpenArea::conf2bot(double w) const {
	auto p = MappedRect::MapToReal(Point(w, 0));
	return std::get<1>(BottomContour().coord_at(p));
}
double RectForOpenArea::top2conf(double w) const {
	auto col = HMCont2D::Contour::WeightPoints(TopContour(), {w});
	return MappedRect::MapToSquare(*col.point(0)).x;
}
double RectForOpenArea::bot2conf(double w) const {
	auto col = HMCont2D::Contour::WeightPoints(BottomContour(), {w});
	return MappedRect::MapToSquare(*col.point(0)).x;
}
double RectForOpenArea::right2conf(double w) const{
	auto col = HMCont2D::Contour::WeightPoints(RightContour(), {w});
	return MappedRect::MapToSquare(*col.point(0)).y;
}
double RectForOpenArea::left2conf(double w) const{
	auto col = HMCont2D::Contour::WeightPoints(LeftContour(), {w});
	return MappedRect::MapToSquare(*col.point(0)).y;
}


//============= CLosed area
RectForClosedArea::RectForClosedArea(const HMCont2D::Contour& side, const HMCont2D::Contour& bottom,
		const HMCont2D::Contour& top):
		MappedRect(side, side, bottom, top, false){
	assert(bottom.is_closed());
	assert(top.is_closed());
	//collect points from top
	vector<Point> ptop; ptop.reserve(top.size());
	for (auto p: top.corner_points()) ptop.push_back(*p);
	if (ptop[0] != *top.first()) ptop.insert(ptop.begin(), *top.first());
	//collect points from bot
	vector<Point> pbot; pbot.reserve(bottom.size());
	for (auto p: bottom.corner_points()) pbot.push_back(*p);
	if (pbot[0] != *bottom.first()) pbot.insert(pbot.begin(), *bottom.first());
	//correct direction
	double a1 = HMCont2D::Contour::Area(top);
	double a2 = HMCont2D::Contour::Area(bottom);
	if (a1<0) std::reverse(ptop.begin(), ptop.end());
	if (a2<0) std::reverse(pbot.begin(), pbot.end());
	//set top, bottom so that top be outer and bot - inner contour
	if (fabs(a1)<fabs(a2)) std::swap(ptop, pbot);
	top_is_outer = (fabs(a1)>fabs(a2));
	//assemble conformal mappint to an annulus
	core = HMMath::Conformal::Annulus::Factory(ptop, pbot);
}

shared_ptr<RectForClosedArea>
RectForClosedArea::Build(const HMCont2D::Contour& bottom, const Point* pstart, double h){
	HMCont2D::Container<HMCont2D::ContourTree> toptree =
		HMCont2D::Contour::Offset(bottom, h, HMCont2D::OffsetTp::CLOSED_POLY);
	//if failed to build single connected area -> do smth
	if (toptree.cont_count() != 1) _THROW_NOT_IMP_;
	auto top = *toptree.nodes[0];
	//find point on a cross between normal from bottom.first() and top contour
	//and set it as a start point for top
	auto first3 = bottom.point_siblings(bottom.first());
	double a = Angle(*first3[2], *first3[1], *first3[0])/2.0;
	Vect q = vecRotate(*first3[0] - *first3[1], a);
	Point np1 = *first3[1], np2 = *first3[1] + q*1000;
	HMCont2D::Contour normal;
	normal.add_value(HMCont2D::Edge(&np1, &np2));
	auto cr = HMCont2D::Contour::Cross(top, normal);
	assert(std::get<0>(cr));
	auto gp = top.GuaranteePoint(std::get<1>(cr), toptree.pdata);
	auto top2 = HMCont2D::Contour::Assemble(top, std::get<1>(gp), std::get<1>(gp));
	//all 4 contours were set. call main routine
	return Build(bottom, top2);
}


shared_ptr<RectForClosedArea>
RectForClosedArea::Build(const HMCont2D::Contour& bottom, const HMCont2D::Contour& top){
	HMCont2D::Contour side;
	side.add_value(HMCont2D::Edge(bottom.first(), top.first()));
	return shared_ptr<RectForClosedArea>(new RectForClosedArea(side, bottom, top));
}



HMCont2D::PCollection RectForClosedArea::MapToReal(const vector<const Point*>& pnt) const{
	//length of curve of inner circle
	double il = 2*M_PI*core->module(); 
	//length of curve of first point of inner circle
	Point p0 = core->InnerCirclePoints()[0];
	double az = atan2(p0.y, p0.x)*core->module();
	//vector of points in canonic area
	vector<Point> cp;
	for (auto p: pnt){
		//l0 - curve length of current point projected to inner circle
		double l0 = p->x*il + az;
		while (l0>il) l0-=il;
		//canonic angle
		double an = l0/(core->module());
		//canonic radius 
		double y = (top_is_outer) ? p->y : 1.0 - p->y;
		double r = y*(1.0 - core->module()) + core->module();
		//point
		cp.push_back(Point(r*cos(an), r*sin(an)));
	}
	//do mapping
	cp = core->MapToOriginal(cp);
	//write answer
	HMCont2D::PCollection ret;
	for (auto p: cp) ret.add_value(p);
	return ret;
}

HMCont2D::PCollection RectForClosedArea::MapToSquare(const vector<const Point*>& pin) const{
	//map to canonic
	vector<Point> pout(pin.size());
	std::transform(pin.begin(), pin.end(), pout.begin(), [](const Point* p){ return *p; });
	pout = core->MapToAnnulus(pout);
	//(x, y) -> (rad, phi)
	std::transform(pout.begin(), pout.end(), pout.begin(), [](const Point& p){
		return Point(sqrt(sqr(p.x) + sqr(p.y)), atan2(p.y, p.x));
	});
	//(rad, phi) -> (bottom length, rad - inner_radius)
	double m = core->module();
	std::transform(pout.begin(), pout.end(), pout.begin(), [&m](const Point& p){
		return Point(p.y*m, p.x - m);
	});
	//normalize
	double maxx = 2*M_PI*m, maxy = 1.0 - m; 
	double x0 = core->module()*core->PhiInner(0);
	std::transform(pout.begin(), pout.end(), pout.begin(), [&](const Point& p){
		Point ret(p.x - x0, p.y);
		if (ret.x < 0) ret.x += 2*M_PI*m;
		ret.x/=maxx; ret.y/=maxy;
		return ret;
	});
	//if top is not out -> reverse y
	if (!top_is_outer){
		std::transform(pout.begin(), pout.end(), pout.begin(), [](const Point& p){
			return Point(p.x, 1.0 - p.y);
		});
	}
	//write answer
	HMCont2D::PCollection ret;
	std::for_each(pout.begin(), pout.end(), [&ret](const Point& p){
		ret.add_value(p);
	});
	return ret;
}

//TODO: move to MappedRect maybe?
double RectForClosedArea::conf2top(double w) const {
	_THROW_NOT_IMP_;
	auto p = MappedRect::MapToReal(Point(w, 1));
	return std::get<1>(TopContour().coord_at(p));
};
double RectForClosedArea::conf2bot(double w) const {
	_THROW_NOT_IMP_;
	auto p = MappedRect::MapToReal(Point(w, 0));
	return std::get<1>(BottomContour().coord_at(p));
}
double RectForClosedArea::top2conf(double w) const {
	_THROW_NOT_IMP_;
	auto col = HMCont2D::Contour::WeightPoints(TopContour(), {w});
	return MappedRect::MapToSquare(*col.point(0)).x;
}
double RectForClosedArea::bot2conf(double w) const {
	auto p = HMCont2D::Contour::WeightPoint(BottomContour(), w);
	return MappedRect::MapToSquare(p).x;
}
double RectForClosedArea::right2conf(double w) const{
	_THROW_NOT_IMP_;
	auto col = HMCont2D::Contour::WeightPoints(RightContour(), {w});
	return MappedRect::MapToSquare(*col.point(0)).y;
}
double RectForClosedArea::left2conf(double w) const{
	_THROW_NOT_IMP_;
	auto col = HMCont2D::Contour::WeightPoints(LeftContour(), {w});
	return MappedRect::MapToSquare(*col.point(0)).y;
}

// ========================== MappedMesher
void MappedMesher::Fill(TBotPart bottom_partitioner, TVertPart vertical_partitioner){
	HMCont2D::Contour bt = rect->BottomContour();
	//1) get real weights from conformal weights at bottom source
	double wbot_start = (ISEQ(wstart, 0.0)) ? 0.0 : rect->conf2bot(wstart);
	double wbot_end = (ISEQ(wend, 1.0)) ? 1.0 : rect->conf2bot(wend);

	//if bottom is one-edge contour with wrong direction
	assert(wbot_start<wbot_end);

	//2) Build bottom partitions
	auto bpart = bottom_partitioner(wbot_start, wbot_end);

	//3) Build partitions of vertical lines
	vector<vector<double>> vlines;
	for (auto v: bpart) vlines.push_back(vertical_partitioner(v));

	int isz = bpart.size();
	int jsz = std::max_element(vlines.begin(), vlines.end(),
		[](vector<double>& x, vector<double>& y){ return x.size() < y.size(); })->size();

	//4) build regular grid.
	GridGeom g4 = GGeom::Constructor::RectGrid01(isz, jsz);
	//5) copy left/right points to pcollections
	vector<int> right_indicies, left_indicies;
	for (int j=0; j<jsz; ++j){
		left_indicies.push_back(j*isz);
		right_indicies.push_back(j*isz + isz - 1);
	}
	left_points.add_values(GGeom::Info::SharePoints(g4, left_indicies));
	right_points.add_values(GGeom::Info::SharePoints(g4, right_indicies));
	//6) delete unused cells
	vector<const Cell*> ucells;
	int kc = 0;
	for (int j=0; j<jsz-1; ++j){
		for (int i=0; i<isz-1; ++i){
			if (vlines[i].size() < j+1 || vlines[i+1].size() < j+1)
				ucells.push_back(g4.get_cell(kc));
			++kc;
		}
	}
	GGeom::Modify::RemoveCells(g4, ucells);
	//7) Remove deleted points from left/right
	left_points.RemoveUnused();
	right_points.RemoveUnused();

	//8) weight coordinates
	for_each(bpart.begin(), bpart.end(), [&](double& x){ x = rect->bot2conf(x); });
	double hx = 1.0/(isz-1);
	double hy = 1.0/(jsz-1);
	auto toweights = [&](GridPoint* p){
		int i = lround(p->x/hx), j = lround(p->y/hy);
		p->set(bpart[i], vlines[i][j]);
	};
	GGeom::Modify::PointModify(g4, toweights);

	//9) modify points using mapping
	vector<const Point*> p(g4.n_points());
	for (int i=0; i<g4.n_points(); ++i) p[i] = g4.get_point(i);
	HMCont2D::PCollection mp = rect->MapToReal(p);
	auto mapfunc = [&mp](GridPoint* p){
		p->x = mp.point(p->get_ind())->x;
		p->y = mp.point(p->get_ind())->y;
	};
	GGeom::Modify::PointModify(g4, mapfunc);

	//10) copy to results
	GGeom::Modify::ShallowAdd(&g4, &result);
}

HMCont2D::Contour MappedMesher::LeftContour(){
	return HMCont2D::Constructor::ContourFromPoints(left_points, false);
}

HMCont2D::Contour MappedMesher::RightContour(){
	return HMCont2D::Constructor::ContourFromPoints(right_points, false);
}








