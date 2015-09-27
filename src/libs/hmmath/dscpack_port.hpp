#ifndef HYBMESH_HMMATH_DSCPACK_PORT_HPP
#define HYBMESH_HMMATH_DSCPACK_PORT_HPP
#include "hmconformal.hpp"

namespace HMMath{ namespace Conformal{ namespace Impl{
namespace DSCPack{

struct Pt{
	double x, y;
};

class ToAnnulus: public HMMath::Conformal::Annulus{
//input. subscripts: 1 for outer, 2 for inner
	const int n1, n2;
	const int prec;
	const vector<Pt> z1, z2;
//output
	vector<double> alfa1, alfa2;
	vector<Pt> w1, w2;
	vector<double> phi1, phi2;
	double u;
	Pt c;
	vector<double> qwork;

//conformal module = rad(inner)
	double _module;

	//prec is number of Gauss-Jacobi integration points.
	//Result will be approximated as O(1E-prec)
	ToAnnulus(const vector<Pt>& outer, const vector<Pt>& inner, int prec);

	//calculates minimum distance between original points in
	//canonic area. Used for detection of the crowding problem.
	double min_wdist() const;
public:
	//all contours are in anti-clockwise direction
	//all points are unique
	//Returns 0 if failed (f.e. due to crowding problem or smth)
	static shared_ptr<ToAnnulus>
	Build(const vector<Point>& outer, const vector<Point>& inner);

	static shared_ptr<ToAnnulus>
	Build(const HMCont2D::Contour& outer, const HMCont2D::Contour& inner);

	double module() const override { return _module; }

	vector<Point> InnerCirclePoints() const override;
	vector<Point> OuterCirclePoints() const override;

	vector<Point> MapToOriginal(const vector<Point>& input) const override;
	vector<Point> MapToAnnulus(const vector<Point>& input) const override;

	double PhiInner(int i) const override;
	double PhiOuter(int i) const override;
};


}
}}}
#endif

