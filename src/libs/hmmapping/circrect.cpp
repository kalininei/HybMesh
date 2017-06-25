#include "circrect.hpp"
#include "rectangle_grid_builder.hpp"
#include "buildcont.hpp"
#include "partcont.hpp"
#include "healgrid.hpp"
#include "buildgrid.hpp"
#include "modgrid.hpp"
#include "assemble2d.hpp"

namespace{

void reflect_and_merge(HM2D::GridData& g, Point p1, Point p2){
	double lx = p2.x - p1.x, ly = p2.y - p1.y;
	double r2 = lx*lx + ly*ly;
	lx /= sqrt(r2); ly /= sqrt(r2);
	double A[4] = {lx*lx - ly*ly, 2*lx*ly, 2*lx*ly, ly*ly - lx*lx};
	HM2D::GridData g2;
	HM2D::DeepCopy(g, g2);
	for (int i=0; i<g2.vvert.size(); ++i){
		auto& p = g2.vvert[i];
		Point r = *p - p1;
		p->x = A[0] * r.x + A[1] * r.y + p1.x;
		p->y = A[2] * r.x + A[3] * r.y + p1.y;
	}
	for (auto& e: g2.vedges){ std::swap(e->left, e->right); }
	for (auto& c: g2.vcells){ std::reverse(c->edges.begin(), c->edges.end()); }
	HM2D::Grid::Algos::MergeTo(g2, g);
}

HM2D::GridData laplace_algo(HM2D::EdgeData& left, HM2D::EdgeData& bot,
		HM2D::EdgeData& right, HM2D::EdgeData& top){
	return HMMap::FDMLaplaceRectGrid(left, bot, right, top);
}

HM2D::GridData linear_algo(HM2D::EdgeData& left, HM2D::EdgeData& bot,
		HM2D::EdgeData& right, HM2D::EdgeData& top){
	HM2D::GridData ret = HM2D::Grid::Constructor::RectGrid01(bot.size(), left.size());

	auto pleft = HM2D::Contour::OrderedPoints(left);
	auto pbot = HM2D::Contour::OrderedPoints(bot);
	auto ptop = HM2D::Contour::OrderedPoints(top);
	auto w = HM2D::Contour::EWeights(left);

	for (int i=0; i<bot.size()+1; ++i){
		Point p1 = *pbot[i];
		Point p2 = *ptop[i];
		for (int j=0; j<left.size()+1; ++j){
			auto gp = ret.vvert[i + pbot.size()*j];
			gp->set(Point::Weigh(p1, p2, w[j]));
		}
	}

	return ret;
}

HM2D::GridData ortho_algo(HM2D::EdgeData& left, HM2D::EdgeData& bot,
		HM2D::EdgeData& right, HM2D::EdgeData& top){
	HM2D::GridData ret = HMMap::OrthogonalRectGrid(left, bot, right, top);
	//modify top points so they lay on the circle
	double rad = HM2D::Contour::Last(left)->y;
	for (int i=0; i<bot.size()+1; ++i){
		int ind = i + left.size()*(bot.size()+1);
		vecSetLen(*ret.vvert[ind], rad);
	}
	return ret;
}

HM2D::GridData ortho_circ_algo(HM2D::EdgeData& left, HM2D::EdgeData& bot,
		HM2D::EdgeData& right, HM2D::EdgeData& top){
	return HMMap::OrthogonalRectGrid(left, top, right, bot);
}
}

HM2D::GridData HMMap::Circ4Prototype(Point center, double rad, int n, std::string algo,
		double a, double hcoef){
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
	auto auxcont = HM2D::Contour::Constructor::FromPoints({0, 0, bav, 0});
	std::map<double, double> auxbasis;
	auxbasis[0] = h1; auxbasis[1] = h2;
	auto auxpart = HM2D::Contour::Algos::WeightedPartition(auxbasis, auxcont);
	vector<double> w(auxpart.size()+1);
	auto wit = w.begin();
	for (auto p: HM2D::Contour::OrderedPoints(auxpart)) *wit++ = p->x/bav;

	//assemble stencil
	//::bottom
	vector<Point> botpoints(n1+1);
	for (int i=0; i<n1+1; ++i) botpoints[i] = Point(a*i/n1, a);
	auto bot = HM2D::Contour::Constructor::FromPoints(botpoints);
	//::top
	vector<Point> toppoints(n1+1);
	for (int i=0; i<n1+1; ++i){
		double alpha = M_PI/2.0 - M_PI/4.0*i/n1;
		toppoints[i] = Point(cos(alpha), sin(alpha));
	}
	auto top = HM2D::Contour::Constructor::FromPoints(toppoints);
	//::left
	Point pleft0(0, a), pleft1(0, 1);
	vector<Point> leftpoints(w.size());
	for (int i=0; i<w.size(); ++i) leftpoints[i] = Point::Weigh(pleft0, pleft1, w[i]);
	auto left = HM2D::Contour::Constructor::FromPoints(leftpoints);
	//::right
	Point pright0(a, a), pright1(sqrt(2)/2.0, sqrt(2)/2.0);
	vector<Point> rightpoints(w.size());
	for (int i=0; i<w.size(); ++i) rightpoints[i] = Point::Weigh(pright0, pright1, w[i]);
	auto right = HM2D::Contour::Constructor::FromPoints(rightpoints);
	
	//make mapping
	HM2D::GridData gcirc = [&](){
		if (algo == "linear") return linear_algo(left, bot, right, top);
		else if (algo == "laplace") return laplace_algo(left, bot, right, top);
		else if (algo == "orthogonal-rect") return ortho_algo(left, bot, right, top);
		else if (algo == "orthogonal-circ") return ortho_circ_algo(left, bot, right, top);
		else throw std::runtime_error("unknown circ4grid algorithm");
	}();

	//orthogonal-circ can give non-uniform rectangle,
	//so we take bottom nodes explicitly from gcirc
	//to use it as a base for internal rectangle grid building
	vector<double> botc(2*n1+1);
	if (algo == "orthogonal-circ"){
		auto botcont = HM2D::Contour::Constructor::CutContour(
			HM2D::Contour::Assembler::GridBoundary1(gcirc),
			Point(0, a), Point(a, a));
		int i=0;
		for (auto p: HM2D::Contour::OrderedPoints(botcont)){
			botc[n1+i] =  p->x;
			botc[n1-i] = -p->x;
			++i;
		}
	} else {
		double h = a/n1;
		for (int i=0; i<2*n1+1; ++i) botc[i] = -a + i*h;
	}

	//enclose grids
	reflect_and_merge(gcirc, pright0, pright1);
	reflect_and_merge(gcirc, pleft0, pleft1);
	reflect_and_merge(gcirc, Point(0, 0), Point(1, 0));

	//inner rectangle grid
	auto grect = HM2D::Grid::Constructor::RectGrid(botc, botc);

	//connect all grids
	HM2D::Grid::Algos::MergeTo(gcirc, grect);

	//scale, translate, return
	for (int i=0; i<grect.vvert.size(); ++i){
		(*grect.vvert[i]) *= rad;
		(*grect.vvert[i]) += center;
	}
	return grect;
}
