#include "dscpack_port.hpp"
#include "contour.hpp"
using namespace HMMap::Conformal::Impl::DSCPack;

extern "C"{

//see comments in dscpack_port.f
double dscpack_init_(
//input
	const int* n1,
	const int* n2,
	const Pt z1[],
	const Pt z2[],
	const int* prec,
//output
	double alfa1[],
	double alfa2[],
	Pt w1[],
	Pt w2[],
	double phi1[],
	double phi2[],
	double *u,
	Pt* c,
	double qwork[]
);

Pt dscpack_backward_(
	const Pt* ww1,
	const int* n1,
	const int* n2,
	const Pt z1[],
	const Pt z2[],
	const int* prec,
	const double alfa1[],
	const double alfa2[],
	const Pt w1[],
	const Pt w2[],
	const double phi1[],
	const double phi2[],
	const double* u,
	const Pt* c,
	const double qwork[]
);

Pt dscpack_forward_(
	const Pt* zz1,
	const int* n1,
	const int* n2,
	const Pt z1[],
	const Pt z2[],
	const int* prec,
	const double alfa1[],
	const double alfa2[],
	const Pt w1[],
	const Pt w2[],
	const double phi1[],
	const double phi2[],
	const double* u,
	const Pt* c,
	const double qwork[]
);

}

ToAnnulus::ToAnnulus(const vector<Pt>& outer, const vector<Pt>& inner, int _prec):
	n1(outer.size()), n2(inner.size()), prec(_prec),
	z1(outer), z2(inner), alfa1(n1), alfa2(n2),
	w1(n1), w2(n2), phi1(n1), phi2(n2),
	qwork(prec*(2*(n1+n2)+3)),
	_module(dscpack_init_(
		&n1, &n2, &z1[0], &z2[0], &prec,
		&alfa1[0], &alfa2[0], &w1[0], &w2[0],
		&phi1[0], &phi2[0], &u, &c, &qwork[0]
	))
{}

shared_ptr<ToAnnulus>
ToAnnulus::Build(const vector<Point>& outer, const vector<Point>& inner){
	//if number of nodes > 40 don't use this method.
	//it is slow and doesn't converge
	if (outer.size() + inner.size() > 40) return 0;

	//building
	vector<Pt> v1(outer.size()), v2(inner.size());
	std::transform(outer.begin(), outer.end(), v1.begin(),
			[](const Point& p){ return Pt{p.x, p.y}; });
	std::transform(inner.begin(), inner.end(), v2.begin(),
			[](const Point& p){ return Pt{p.x, p.y}; });

	shared_ptr<ToAnnulus> ret(new ToAnnulus(v1, v2, 12));
	//if points in canonic area are distinguishable everything is ok
	if (ret->min_wdist() > 100*geps) return ret;

	//failed to build -> return nothing
	return 0;
}

shared_ptr<ToAnnulus>
ToAnnulus::Build(const HM2D::EdgeData& outer, const HM2D::EdgeData& inner){
	assert(HM2D::Contour::IsClosed(outer) && HM2D::Contour::IsClosed(inner));
	vector<Point> v1(outer.size()), v2(inner.size());
	auto p1 = HM2D::Contour::OrderedPoints(outer); p1.pop_back();
	auto p2 = HM2D::Contour::OrderedPoints(inner); p2.pop_back();
	std::transform(p1.begin(), p1.end(), v1.begin(), [](shared_ptr<HM2D::Vertex> p){
			return *p;});
	std::transform(p2.begin(), p2.end(), v2.begin(), [](shared_ptr<HM2D::Vertex> p){
			return *p;});
	return Build(v1, v2);
}

double ToAnnulus::min_wdist() const{
	double ret = 1e10;
	//outer polygon
	for (int i=0; i<n1; ++i){
		int inext = (i==n1-1)?0:i+1;
		double d2 = sqr(w1[i].x-w1[inext].x) + sqr(w1[i].y-w1[inext].y);
		if (d2<ret) ret = d2;
	}
	//inner polygon
	for (int i=0; i<n2; ++i){
		int inext = (i==n2-1)?0:i+1;
		double d2 = sqr(w2[i].x-w2[inext].x) + sqr(w2[i].y-w2[inext].y);
		if (d2<ret) ret = d2;
	}
	return sqrt(ret);
}

vector<Point> ToAnnulus::MapToOriginal(const vector<Point>& input) const{
	vector<Point> ret; ret.reserve(input.size());
	auto add_point = [&ret](const Pt& p){ ret.push_back(Point(p.x, p.y)); };

	for (auto& p: input){
		Pt x {p.x, p.y};
		//search amoung originals
		auto fnd1 = std::find_if(w1.begin(), w1.end(), [&x](const Pt& wp){
					return ISEQ(x.x, wp.x) && ISEQ(x.y, wp.y);}); 
		if (fnd1 != w1.end()) {
			add_point( z1[fnd1 - w1.begin()] );
			continue;
		}
		auto fnd2 = std::find_if(w2.begin(), w2.end(), [&x](const Pt& wp){
					return ISEQ(x.x, wp.x) && ISEQ(x.y, wp.y);}); 
		if (fnd2 != w2.end()) {
			add_point( z2[fnd2 - w2.begin()] );
			continue;
		}
		//do mapping
		x = dscpack_backward_(
			&x,
			&n1,
			&n2,
			&z1[0],
			&z2[0],
			&prec,
			&alfa1[0],
			&alfa2[0],
			&w1[0],
			&w2[0],
			&phi1[0],
			&phi2[0],
			&u,
			&c,
			&qwork[0]
		);
		add_point(x);
	}
	return ret;
}

vector<Point> ToAnnulus::MapToAnnulus(const vector<Point>& input) const{
	vector<Point> ret; ret.reserve(input.size());
	auto add_point = [&ret](const Pt& p){ ret.push_back(Point(p.x, p.y)); };

	for (auto& p: input){
		Pt x {p.x, p.y};
		//search amoung originals
		auto fnd1 = std::find_if(z1.begin(), z1.end(), [&x](const Pt& wp){
			return ISEQ(x.x, wp.x) && ISEQ(x.y, wp.y);}); 
		if (fnd1 != z1.end()) {
			add_point( w1[fnd1 - z1.begin()] );
			continue;
		}
		auto fnd2 = std::find_if(z2.begin(), z2.end(), [&x](const Pt& wp){
			return ISEQ(x.x, wp.x) && ISEQ(x.y, wp.y);}); 
		if (fnd2 != z2.end()) {
			add_point( w2[fnd2 - z2.begin()] );
			continue;
		}
		//do mapping
		x = dscpack_forward_(
			&x,
			&n1,
			&n2,
			&z1[0],
			&z2[0],
			&prec,
			&alfa1[0],
			&alfa2[0],
			&w1[0],
			&w2[0],
			&phi1[0],
			&phi2[0],
			&u,
			&c,
			&qwork[0]
		);
		add_point(x);
	}
	return ret;
}

vector<Point> ToAnnulus::InnerCirclePoints() const{
	vector<Point> ret(n2);
	std::transform(w2.begin(), w2.end(), ret.begin(),
		[](const Pt& p){ return Point(p.x, p.y); });
	return ret;
}

vector<Point> ToAnnulus::OuterCirclePoints() const{
	vector<Point> ret(n1);
	std::transform(w1.begin(), w1.end(), ret.begin(),
		[](const Pt& p){ return Point(p.x, p.y); });
	return ret;
}

double ToAnnulus::PhiInner(int i) const{ return phi2[i]; }

double ToAnnulus::PhiOuter(int i) const{ return phi1[i]; }

