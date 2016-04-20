#include "circrect.hpp"
#include "rectangle_grid_builder.hpp"
#include "procgrid.h"

namespace{

void reflect_and_merge(GridGeom& g, Point p1, Point p2){
	double lx = p2.x - p1.x, ly = p2.y - p1.y;
	double r2 = lx*lx + ly*ly;
	lx /= sqrt(r2); ly /= sqrt(r2);
	double A[4] = {lx*lx - ly*ly, 2*lx*ly, 2*lx*ly, ly*ly - lx*lx};
	auto reflect_point=[&p1, &A](GridPoint* p){
		Point r = *p - p1;
		p->x = A[0] * r.x + A[1] * r.y + p1.x;
		p->y = A[2] * r.x + A[3] * r.y + p1.y;
	};
	auto g2 = GGeom::Constructor::DeepCopy(g);
	GGeom::Modify::PointModify(g2, reflect_point);
	GGeom::Modify::ShallowAdd(&g2, &g);
	GGeom::Repair::Heal(g);
}

};

GridGeom HMGMap::Circ4Prototype(Point center, double rad, int n, double a, double hcoef){
	if (a>1.4) throw std::runtime_error("square side is too big");

	assert(n % 8 == 0);
	int n1 = n/8; n = n1*8;
	a/=2;

	//lengths of rays
	double b1 = 1.0 - a;
	double b2 = 1.0 - sqrt(2.0)/2.0;

	//steps at square, outer circle boundary
	double h1 = a/n1;
	double h2 = hcoef * M_PI/4.0/n1;

	//weights for left, right contours
	//double bav = (b1 + b2)/2.0;
	double bav = b1;
	auto auxcont = HMCont2D::Constructor::ContourFromPoints({0, 0, bav, 0});
	std::map<double, double> auxbasis;
	auxbasis[0] = h1; auxbasis[1] = h2;
	auto auxpart = HMCont2D::Algos::WeightedPartition(auxbasis, auxcont, auxcont.pdata);
	vector<double> w(auxpart.size()+1);
	auto wit = w.begin();
	for (auto p: auxpart.ordered_points()) *wit++ = p->x/bav;

	//assemble stencil
	//::bottom
	vector<Point> botpoints(n1+1);
	for (int i=0; i<n1+1; ++i) botpoints[i] = Point(a*i/n1, a);
	auto bot = HMCont2D::Constructor::ContourFromPoints(botpoints);
	//::top
	vector<Point> toppoints(n1+1);
	for (int i=0; i<n1+1; ++i){
		double alpha = M_PI/2.0 - M_PI/4.0*i/n1;
		toppoints[i] = Point(cos(alpha), sin(alpha));
	}
	auto top = HMCont2D::Constructor::ContourFromPoints(toppoints);
	//::left
	Point pleft0(0, a), pleft1(0, 1);
	vector<Point> leftpoints(w.size());
	for (int i=0; i<w.size(); ++i) leftpoints[i] = Point::Weigh(pleft0, pleft1, w[i]);
	auto left = HMCont2D::Constructor::ContourFromPoints(leftpoints);
	//::right
	Point pright0(a, a), pright1(sqrt(2)/2.0, sqrt(2)/2.0);
	vector<Point> rightpoints(w.size());
	for (int i=0; i<w.size(); ++i) rightpoints[i] = Point::Weigh(pright0, pright1, w[i]);
	auto right = HMCont2D::Constructor::ContourFromPoints(rightpoints);
	
	//make mapping
	auto gcirc = HMGMap::FDMLaplasRectGrid(left, bot, right, top);

	//enclose grids
	reflect_and_merge(gcirc, pright0, pright1);
	reflect_and_merge(gcirc, pleft0, pleft1);
	reflect_and_merge(gcirc, Point(0, 0), Point(1, 0));

	//inner rectangle grid
	auto grect = GGeom::Constructor::RectGrid(
			Point(-a, -a), Point(a, a), 2*n1, 2*n1);

	//connect all grids
	GGeom::Modify::ShallowAdd(&gcirc, &grect);
	GGeom::Repair::Heal(grect);

	//scale, translate, return
	auto pmod = [&rad, &center](GridPoint* p){
		p->x *= rad; p->x += center.x;
		p->y *= rad; p->y += center.y;
	};
	GGeom::Modify::PointModify(grect, pmod);
	return grect;
}
