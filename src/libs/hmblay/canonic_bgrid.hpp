#ifndef HYBMESH_HMBLAY_CANONIC_HPP
#define HYBMESH_HMBLAY_CANONIC_HPP
#include "bgrid.hpp"
#include "hmconformal.hpp"

namespace HMBlay{ namespace Impl{

//takes two closed contours (c1, c2) and builds
//averaged one which goes from p1 (in c1) to p2 (in c2)
HM2D::EdgeData
ContoursWeight(const HM2D::EdgeData& c1, Point p1,
		const HM2D::EdgeData& c2, Point p2);


class MappedRect{
	bool _use_rect_approx;
protected:
	HM2D::EdgeData left, right, bottom, top;

	MappedRect(const HM2D::EdgeData& _left, const HM2D::EdgeData& _right,
			const HM2D::EdgeData& _bottom, const HM2D::EdgeData& _top,
			bool rect_approx): _use_rect_approx(rect_approx){
		HM2D::DeepCopy(_left, left);
		HM2D::DeepCopy(_right, right);
		HM2D::DeepCopy(_bottom, bottom);
		HM2D::DeepCopy(_top, top);
	}
public:
	bool use_rect_approx() const { return _use_rect_approx; }
	//real contour mapped by {u=0->1, v=1 }
	HM2D::EdgeData TopContour() const { return HM2D::EdgeData(top); }
	//real contour mapped by {u=0->1, v=0 }
	HM2D::EdgeData BottomContour() const { return HM2D::EdgeData(bottom); }
	//real contour mapped by {u=0, v=0->1 }
	HM2D::EdgeData LeftContour() const { return HM2D::EdgeData(left); }
	//real contour mapped by {u=1, v=0->1 }
	HM2D::EdgeData RightContour() const { return HM2D::EdgeData(right); }

	//real coordinates from normalized in [0,1] square
	virtual HM2D::VertexData MapToReal(const vector<const Point*>& p) const = 0;
	HM2D::VertexData MapToReal(const vector<Point*>& p) const{
		return MapToReal(vector<const Point*>(p.begin(), p.end()));
	}
	Point MapToReal(const Point& p) const{ return *(MapToReal({&p})[0]); }
	
	//square coordinates from real coordinates
	virtual HM2D::VertexData MapToSquare(const vector<const Point*>& p) const = 0;
	HM2D::VertexData MapToSquare(const vector<Point*>& p) const{
		return MapToSquare(vector<const Point*>(p.begin(), p.end()));
	}
	Point MapToSquare(const Point& p) const { return *(MapToSquare({&p})[0]); }

	virtual Point MapBndToReal(const Point& p) const { return MapToReal(p); }
	virtual Point MapBndToSquare(const Point& p) const { return MapToSquare(p); }
	
	//weights transformers
	virtual double conf2top(double w) const = 0;
	virtual double conf2bot(double w) const = 0;
	virtual double top2conf(double w) const = 0;
	virtual double bot2conf(double w) const = 0;
	virtual double right2conf(double w) const = 0;
	virtual double left2conf(double w) const = 0;
	double top2bot(double w) const;
	double bot2top(double w) const;

	// ================ constructors:
	//1. From four lines
	//direction of contours:
	//left   (u=0;    v: 0->1)
	//right  (u=1;    v: 0->1)
	//bottom (u=0->1; v: 0)
	//top    (u=0->1; v: 1)
	static shared_ptr<MappedRect>
	Factory(HM2D::EdgeData& left, HM2D::EdgeData& right,
			HM2D::EdgeData& bottom, HM2D::EdgeData& top,
			int femn, bool use_rect_approx);
	//build forcing usage of rectangle approximation
	static shared_ptr<MappedRect>
	RectApprox(HM2D::EdgeData& left, HM2D::EdgeData& right,
			HM2D::EdgeData& bottom, HM2D::EdgeData& top);

	//2. From three lines (no top) and vertical distance h
	static shared_ptr<MappedRect>
	Factory(HM2D::EdgeData& left, HM2D::EdgeData& right,
			HM2D::EdgeData& bottom, double h,
			int femn, bool use_rect_approx);
};

class RectForOpenArea: public MappedRect{
protected:
	shared_ptr<HMMap::Conformal::Rect> core;
public:
	RectForOpenArea(HM2D::EdgeData& left, HM2D::EdgeData& right,
			HM2D::EdgeData& bottom, HM2D::EdgeData& top,
			int femn, bool use_rect_approx, bool force_rect_approx=false);
	HM2D::VertexData MapToReal(const vector<const Point*>& p) const override;
	HM2D::VertexData MapToSquare(const vector<const Point*>& p) const override;

	double conf2top(double w) const override;
	double conf2bot(double w) const override;
	double top2conf(double w) const override;
	double bot2conf(double w) const override;
	double right2conf(double w) const override;
	double left2conf(double w) const override;
};

class RectForClosedArea: public MappedRect{
	shared_ptr<HMMap::Conformal::Annulus> core;
	bool top_is_outer;
public:
	//if |area(bottom)|<|area(top)| => contours have negative direction
	//else => positive
	RectForClosedArea(const HM2D::EdgeData& side, const HM2D::EdgeData& bottom,
			const HM2D::EdgeData& top);

	//pstart should lie in bottom
	//if bottom is in anti-clockwise direction -> top is an outer contour
	static shared_ptr<RectForClosedArea>
	Build(const HM2D::EdgeData& bottom, const Point* pstart, double h);

	//bottom.first(), top.first() will be connected
	static shared_ptr<RectForClosedArea>
	Build(const HM2D::EdgeData& bottom, const HM2D::EdgeData& top);

	//mapping functions
	HM2D::VertexData MapToReal(const vector<const Point*>& p) const override;
	HM2D::VertexData MapToSquare(const vector<const Point*>& p) const override;

	Point MapBndToReal(const Point& p) const override;
	Point MapBndToSquare(const Point& p) const override;

	double conf2top(double w) const override;
	double conf2bot(double w) const override;
	double top2conf(double w) const override;
	double bot2conf(double w) const override;
	double right2conf(double w) const override;
	double left2conf(double w) const override;
};


// ============================== Meshers
class MappedMesher{
public:
	//start and end square x coordinates of meshing
	double wstart, wend;
	MappedRect* rect;
	MappedMesher(MappedRect* r, double w1, double w2){
		rect = r; wstart = w1; wend = w2;
	};

	vector<GridPoint*> right_points, left_points;
	BGrid result;

	//takes weight start, weight end of source line.
	//returns vector of desired partition points by weight coordinate
	//including start and end weights
	typedef std::function<vector<double>(double, double)> TBotPart;
	//takes point by weight coordinate on source line,
	//return weight depths of blayers starting from 0.0
	//depth = 1 is the height of rect at certain weight
	typedef std::function<vector<double>(double)>  TVertPart;

	//Build a grid.
	//source is used for weight calculation
	// 1 - left, 2 - bottom, 3 - right, 4 - top
	void Fill(TBotPart bottom_partitioner, TVertPart vertical_partitioner, int source);

	//returns grid boundaries as shallow copy
	HM2D::EdgeData LeftContour();
	HM2D::EdgeData RightContour();
	
};



}}
#endif

