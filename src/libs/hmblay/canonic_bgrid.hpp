#ifndef HYBMESH_HMBLAY_CANONIC_HPP
#define HYBMESH_HMBLAY_CANONIC_HPP
#include "bgrid.hpp"

namespace HMBlay{ namespace Impl{


class MappedRect{
public:
	//real contour mapped by {u=0->1, v=1 }
	virtual HMCont2D::Contour TopContour() = 0;
	//real contour mapped by {u=0->1, v=0 }
	virtual HMCont2D::Contour BottomContour() = 0;
	//real contour mapped by {0,0 -> 1,0 -> 1,1 -> 0,1}
	virtual HMCont2D::Contour FullContour() = 0;
	virtual HMCont2D::PCollection MapToReal(const vector<const Point*>& p) = 0;
};

// ================= FourLineRect: all 4 lines are known
class FourLineRect: public MappedRect{
public:
	//direction of contours:
	//left   (u=0;    v: 0->1)
	//right  (u=1;    v: 0->1)
	//bottom (u=0->1; v: 0)
	//top    (u=0->1; v: 1)
	static shared_ptr<FourLineRect>
	Factory(HMCont2D::Contour& left, HMCont2D::Contour& right,
			HMCont2D::Contour& bottom, HMCont2D::Contour& top);
};

// ================== Cavity: unknown upper bound
class MappedCavity: public MappedRect{
public:
	static shared_ptr<MappedCavity>
	Factory(HMCont2D::Contour& left, HMCont2D::Contour& right,
			HMCont2D::Contour& bottom, double h);
};


class FourStraightLinesRect: public MappedCavity, public FourLineRect{
	std::array<Point, 4> pts;
public:
	//topleft, bottomleft, bottomright, topright
	FourStraightLinesRect(Point p0, Point p1, Point p2, Point p3): pts{p0, p1, p2, p3} {};

	//Overridden
	HMCont2D::Contour TopContour() override;
	HMCont2D::Contour BottomContour() override;
	HMCont2D::Contour FullContour() override;
	HMCont2D::PCollection MapToReal(const vector<const Point*>& p) override;
};

// ============================== Meshers
class MappedMesher{
public:
	double wstart, wend;
	MappedRect* rect;
	MappedMesher(MappedRect* r, double w1, double w2){
		rect = r; wstart = w1; wend = w2;
	};

	HMCont2D::PCollection right_points, left_points;
	BGrid result;

	//takes length start, length end of source line.
	//returns vector of desired partition points by length coordinate
	//including length start and length end
	typedef std::function<vector<double>(double, double)> TBotPart;
	//takes point by length coordinate on source line,
	//return real depths of blayers starting from 0.0
	typedef std::function<vector<double>(double)>  TVertPart;

	//Build a grid
	void Fill(TBotPart bottom_partitioner, TVertPart vertical_partitioner);

	//returns grid boundaries as shallow copy
	HMCont2D::Contour LeftContour();
	HMCont2D::Contour RightContour();
};



}}
#endif

