#include "canonic_bgrid.hpp"
#include "hmfem.hpp"

using namespace HMBlay::Impl;

shared_ptr<MappedCavity>
MappedCavity::Factory(HMCont2D::Contour& left, HMCont2D::Contour& right,
		HMCont2D::Contour& bottom, double h){
	bool straight_left = ISEQ(left.length(), Point::dist(*left.first(), *left.last()));
	bool straight_right = ISEQ(right.length(), Point::dist(*right.first(), *right.last()));
	bool straight_bottom = ISEQ(bottom.length(), Point::dist(*bottom.first(), *bottom.last()));
	//if all boundaries are straight lines
	if (straight_left && straight_right && straight_bottom){
		Point p1 = *left.first();
		Point p3 = *right.first();
		Vect v1 = *left.last() - *left.first(); vecNormalize(v1);
		Vect v2 = *right.last() - *right.first(); vecNormalize(v2);
		return std::make_shared<FourStraightLinesRect>(p1 + v1*h, p1, p3, p3 + v2*h);
	}
	_THROW_NOT_IMP_;
}

shared_ptr<FourLineRect>
FourLineRect::Factory(HMCont2D::Contour& left, HMCont2D::Contour& right,
		HMCont2D::Contour& bottom, HMCont2D::Contour& top){
	assert(*left.first() == *bottom.first());
	assert(*left.last() == *top.first());
	assert(*right.first() == *bottom.last());
	assert(*right.last() == *top.last());

	bool straight_left = ISEQ(left.length(), Point::dist(*left.first(), *left.last()));
	bool straight_right = ISEQ(right.length(), Point::dist(*right.first(), *right.last()));
	bool straight_bottom = ISEQ(bottom.length(), Point::dist(*bottom.first(), *bottom.last()));
	bool straight_top = ISEQ(top.length(), Point::dist(*top.first(), *top.last()));
	if (straight_left && straight_right && straight_top && straight_bottom){
		Point p1 = *left.last();
		Point p2 = *left.first();
		Point p3 = *right.first();
		Point p4 = *right.last();
		return std::make_shared<FourStraightLinesRect>(p1, p2, p3, p4);
	}
	_THROW_NOT_IMP_;
};


HMCont2D::Contour FourStraightLinesRect::TopContour() {
	HMCont2D::Contour ret;
	ret.add_value(HMCont2D::Edge(&pts[0], &pts[3]));
	return ret;
}

HMCont2D::Contour FourStraightLinesRect::BottomContour() {
	HMCont2D::Contour ret;
	ret.add_value(HMCont2D::Edge(&pts[1], &pts[2]));
	return ret;
}
HMCont2D::Contour FourStraightLinesRect::FullContour() {
	HMCont2D::Contour ret;
	ret.add_value(HMCont2D::Edge(&pts[1], &pts[2]));
	ret.add_value(HMCont2D::Edge(&pts[2], &pts[3]));
	ret.add_value(HMCont2D::Edge(&pts[3], &pts[0]));
	ret.add_value(HMCont2D::Edge(&pts[0], &pts[1]));
	return ret;
}

HMCont2D::PCollection FourStraightLinesRect::MapToReal(const vector<const Point*>& points){
	HMCont2D::PCollection ret;
	for (auto p: points){
		Point pup   = Point::Weigh(pts[0], pts[3], p->x);
		Point pdown = Point::Weigh(pts[1], pts[2], p->x);
		ret.add_value(Point::Weigh(pdown, pup, p->y));
	}
	return ret;
}

void MappedMesher::Fill(TBotPart bottom_partitioner, TVertPart vertical_partitioner){
	//0) Build partitions
	HMCont2D::Contour bt = rect->BottomContour();
	double blen = bt.length();
	auto bpart = bottom_partitioner(wstart*blen, wend*blen);
	vector<vector<double>> vlines;
	for (auto v: bpart) vlines.push_back(vertical_partitioner(v));
	int isz = bpart.size();
	int jsz = std::max_element(vlines.begin(), vlines.end(),
		[](vector<double>& x, vector<double>& y){ return x.size() < y.size(); })->size();
	double jmax = std::max_element(vlines.begin(), vlines.end(),
		[](vector<double>& x, vector<double>& y){ return x.back() < y.back(); })->back();

	//1) build regular grid.
	GridGeom g4 = GGeom::Constructor::RectGrid01(isz, jsz);
	//2) copy left/right points to pcollections
	vector<int> right_indicies, left_indicies;
	for (int j=0; j<jsz; ++j){
		left_indicies.push_back(j*isz);
		right_indicies.push_back(j*isz + isz - 1);
	}
	left_points.add_values(GGeom::Info::SharePoints(g4, left_indicies));
	right_points.add_values(GGeom::Info::SharePoints(g4, right_indicies));

	//2) delete unused cells
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

	//4) Remove deleted points from left/right
	left_points.RemoveUnused();
	right_points.RemoveUnused();

	//3) stretch grid points according to partition and fit to [0, 1] square
	auto stretchfun = [&](GridPoint* p){
		int i = std::lround(p->x*(isz-1));
		int j = std::lround(p->y*(jsz-1));
		p->x = bpart[i];
		p->y = vlines[i][j];
		p->x/=blen;
		p->y/=jmax;
	};
	GGeom::Modify::PointModify(g4, stretchfun);
	
	//5) modify points using mapping
	vector<const Point*> p(g4.n_points());
	for (int i=0; i<g4.n_points(); ++i) p[i] = g4.get_point(i);
	HMCont2D::PCollection mp = rect->MapToReal(p);
	auto mapfunc = [&mp](GridPoint* p){
		p->x = mp.point(p->get_ind())->x;
		p->y = mp.point(p->get_ind())->y;
	};
	GGeom::Modify::PointModify(g4, mapfunc);

	//6) copy to results
	GGeom::Modify::ShallowAdd(&g4, &result);
}

HMCont2D::Contour MappedMesher::LeftContour(){
	return HMCont2D::Constructor::ContourFromPoints(left_points, false);
}

HMCont2D::Contour MappedMesher::RightContour(){
	return HMCont2D::Constructor::ContourFromPoints(right_points, false);
}








