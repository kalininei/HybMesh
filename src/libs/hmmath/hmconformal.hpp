#ifndef HYBMESH_HMMATH_HMCONFORMAL_HPP
#define HYBMESH_HMMATH_HMCONFORMAL_HPP
#include "hmproject.h"
#include "bgeom2d.h"
#include "hybmesh_contours2d.hpp"
#include "hmcallback.hpp"

namespace HMMath{ namespace Conformal{

struct Options{
	//Options are used for both Rect and Annulus constructors
	Options():
		use_rect_approx(true),     //rectangle approximation possibility
		use_scpack(true),          //scpack/dscpack solution possibility (for rect/annulus)
		right_angle_eps(geps),     //negligible right angle diviation
		length_weight(1.0+geps),   //negligible lenght deviation at which polyline is treated as a straight line
		scpack_ratio_limit(5.1),   //width/height limit at which scpack is used
		fem_nrec(1000){}           //recommended number of fem points

	// ===== rectangle approximation settings
	//if we consider rectangle approximation
	bool use_rect_approx;
	//semi-analytical solution possibility (solvers are very unstable and slow)
	bool use_scpack;
	//if side angles are outside (pi/2 +/- right_angle_eps) => ignore rect approx
	double right_angle_eps;
	//if line.length() > length_weight * (end-start) => ignore rect approx
	double length_weight;
	//estimate the possibility of using rect approx for given polygon
	bool CanUseRectApprox(const vector<Point>& path, int i1, int i2, int i3) const;

	// ===== scpack setting
	//if module > scpack_ratio_limit => scpack is ignored due to crowding problem
	double scpack_ratio_limit;
	//estimate the possibility of using scpack for given polygon
	bool CanUseSCPACK(const vector<Point>& path, int i1, int i2, int i3) const;

	// ===== dcpack setting
	//if module > scpack_ratio_limit => scpack is ignored due to crowding problem
	bool CanUseDSCPACK() const;

	//===== fem settings
	//recommended fem grid partition
	int fem_nrec;
};


//Mapping from closed polygon to rectangle [0, h]x[0,1]
//h - conformal module of input polygon
class Rect{
public:
	virtual double module() const = 0;
	virtual vector<Point> MapToPolygon(const vector<Point>& input) const = 0;
	virtual vector<Point> MapToRectangle(const vector<Point>& input) const = 0;
	virtual vector<Point> MapToPolygon(const vector<Point*>& input) const;
	virtual vector<Point> MapToRectangle(const vector<Point*>& input) const;

	Point MapToPolygon1(Point p) const { return MapToPolygon(vector<Point> {p})[0]; }
	Point MapToRectangle1(Point p) const { return MapToRectangle(vector<Point> {p})[0]; }

	//get canonic rectangle containing all original points
	virtual vector<Point> RectPoints() const = 0;
	HMCont2D::Container<HMCont2D::Contour> RectContour() const;

	//path -- an ordered set of points in anti-clockwise direction
	//corners -- indicies or corner points directed as
	//   top-left -> bottom-left -> bottom-right -> bottom-top
	static shared_ptr<Rect> Factory(
			const vector<Point>& path,
			std::array<int, 4> corners,
			const Options& opt = Options());

	static shared_ptr<Rect> Factory(
			const HMCont2D::Contour& left,
			const HMCont2D::Contour& right,
			const HMCont2D::Contour& bot,
			const HMCont2D::Contour& top,
			const Options& opt = Options());

	//prepare data for factory input from 4 connected contours
	//directed in such a way that
	//area is bounded in an anti-clockwise direction by
	//  reverse(left), bottom, right, reverse(top)
	//  returning std::array starts with zero
	static std::tuple<vector<Point>, std::array<int, 4>>
	FactoryInput(const HMCont2D::Contour& left, const HMCont2D::Contour& right,
		const HMCont2D::Contour& bot, const HMCont2D::Contour& top);

	//from closed contour
	static std::tuple<vector<Point>, std::array<int, 4>>
	FactoryInput(const HMCont2D::Contour& cont, const std::array<int,4>& a);
};

struct TBuildRect: public HMCallback::ExecutorBase{
	HMCB_SET_PROCNAME("Conformal mapping to rectangle");
	HMCB_SET_DEFAULT_DURATION(100);

	shared_ptr<Rect> _run(const HMCont2D::Contour& left,
		const HMCont2D::Contour& right,
		const HMCont2D::Contour& bot,
		const HMCont2D::Contour& top,
		const Options& opt);
};
extern HMCallback::FunctionWithCallback<TBuildRect> BuildRect;

class Annulus{
public:
	//points should be ordered in anti-clockwise direction
	//inner should lie strictly within outer
	//all points should be unique
	static shared_ptr<Annulus> Factory(
		const vector<Point>& outerpnt,
		const vector<Point>& innerpnt,
		const Options& opt = Options());

	//= RadInner < 1.0
	virtual double module() const = 0;
	//get original points mapped to inner circle
	virtual vector<Point> InnerCirclePoints() const = 0;
	//get original points mapped to outer circle
	virtual vector<Point> OuterCirclePoints() const = 0;
	//mapping functions
	virtual vector<Point> MapToOriginal(const vector<Point>& input) const = 0;
	virtual vector<Point> MapToAnnulus(const vector<Point>& input) const = 0;


	Point MapToOriginal1(const Point& p) const { return MapToOriginal(vector<Point> {p})[0]; }
	Point MapToAnnulus1(const Point& p) const { return MapToAnnulus(vector<Point> {p})[0]; }
	
	virtual Point MapToAnnulusBnd(const Point& p) const { return MapToAnnulus1(p); }
	virtual Point MapToOriginalBnd(const Point& p) const { return MapToOriginal1(p); }

	//calculate angle at which i-th point is mapped to annulus
	virtual double PhiInner(int i) const = 0;
	virtual double PhiOuter(int i) const = 0;
	
	HMCont2D::Container<HMCont2D::Contour> InnerCircleContour() const;
	HMCont2D::Container<HMCont2D::Contour> OuterCircleContour() const;
};


namespace Impl{

//Approximatin of conformal mapping to rectangle using simple geometry
class RectApprox: public Rect{
	double _len_top, _len_bot, _len_left, _len_right;
	double _module;
	bool _left_straight, _right_straight, _top_straight, _bot_straight;
	HMCont2D::PCollection pcol;
	HMCont2D::Contour left, right, top, bot;

	//which direction will be used for approximation
	//x-ksi if true and y-eta otherwise
	bool _xksi;
	
	vector<Point> MapToPolygon_xksi(const vector<Point>& input) const;
	vector<Point> MapToPolygon_yeta(const vector<Point>& input) const;
	Point MapToRectangle_xksi(const Point& input) const;
	Point MapToRectangle_yeta(const Point& input) const;
public:

	double module() const override { return _module; }
	vector<Point> MapToPolygon(const vector<Point>& input) const override;
	Point MapToPolygon(const Point& input) const;

	//mapping functions
	vector<Point> MapToRectangle(const vector<Point>& input) const override;
	Point MapToRectangle(const Point& p) const;

	//Builders
	static shared_ptr<RectApprox> Build(const vector<Point>& path, int i1, int i2, int i3);
	static shared_ptr<RectApprox>
	Build(const HMCont2D::Contour& left, const HMCont2D::Contour& right,
		const HMCont2D::Contour& bot, const HMCont2D::Contour& top);

	//get mapping of original points
	vector<Point> RectPoints() const override {_THROW_NOT_IMP_;}
};

}
}}
#endif
