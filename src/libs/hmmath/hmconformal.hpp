#ifndef HYBMESH_HMMATH_HMCONFORMAL_HPP
#define HYBMESH_HMMATH_HMCONFORMAL_HPP
#include "hmproject.h"
#include "bgeom2d.h"
#include "hybmesh_contours2d.hpp"

namespace HMMath{ namespace Conformal{

//Mapping from closed polygon to rectangle [0, h]x[0,1]
//h - conformal module of input polygon
class Rect{
public:
	struct Option{
		//if module > SCPACK_RATIO_LIMIT => scpack is ignored
		//due to crowding problem
		const double SCPACK_RATIO_LIMIT = 5.1;

		Option():
			use_rect_approx(true),
			right_angle_eps(geps),
			length_weight(1.0+geps){}
		//rectangle approximation settings
		bool use_rect_approx;
		double right_angle_eps;
		double length_weight;
		bool CanUseRectApprox(const vector<Point>& path, int i1, int i2, int i3) const;
		//elongate strip areas ddm approximation settings
		//...
		bool CanUseSCPACK(const vector<Point>& path, int i1, int i2, int i3) const;
	};

	virtual double module() const = 0;
	virtual vector<Point> MapToPolygon(const vector<Point>& input) const = 0;
	virtual vector<Point> MapToRectangle(const vector<Point>& input) const =0;
	//get rectangle
	virtual vector<Point> RectPoints() const = 0;
	HMCont2D::Container<HMCont2D::Contour> RectContour() const;

	//path -- an ordered set of points in anti-clockwise direction
	//corners -- indicies or corner points directed as
	//   top-left -> bottom-left -> bottom-right -> bottom-top
	static shared_ptr<Rect> Factory(
			const vector<Point>& path,
			std::array<int, 4> corners,
			const Option& opt = Option());

	//prepare data for factory input from 4 connected contours
	//directed in such a way that
	//area is bounded in an anti-clockwise direction by
	//  reverse(left), bottom, right, reverse(top)
	//  returning std::array starts with zero
	static std::tuple<vector<Point>, std::array<int, 4>>
	FactoryInput(const HMCont2D::Contour& left, const HMCont2D::Contour& right,
		const HMCont2D::Contour& bot, const HMCont2D::Contour& top);
};


namespace Impl{

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
