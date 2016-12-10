#ifndef HYBMESH_BGEOM_H
#define HYBMESH_BGEOM_H

#include "hmproject.h"
#include <algorithm>
#include "addalgo.hpp"

//Point
struct Point{
	double x,y;
	//constructors
	explicit Point(double _x=0, double _y=0) noexcept:x(_x), y(_y){}
	Point(const Point& p2) noexcept:x(p2.x), y(p2.y){}
	Point& operator=(const Point& p) noexcept{
		if (this!=&p){
			this->x=p.x; this->y=p.y;
		}
		return *this;
	}
	void set(double a=0, double b=0) noexcept {x=a; y=b;}
	void set(const Point& p2) noexcept {x=p2.x; y=p2.y;}
	//from plain arrays
	template<int Dim>
	static vector<Point> read_from_plain(const vector<double>& pcoord){
		vector<Point> ret; ret.reserve(pcoord.size()/Dim);
		auto it=pcoord.begin();
		while (it!=pcoord.end()){
			double x= *it++;
			double y=(Dim>1) ? *it++ : 0.0;
			ret.push_back(Point(x,y));
		}
		return ret;
	}
	// ============== measures and distances
	static double meas(const Point& p1, const Point& p2) noexcept{
		return (p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y);
	}
	static double dist(const Point& p1, const Point& p2) noexcept{
		return sqrt(meas(p1,p2));
	}
	//measure from p to [L1, L2]. ksi is the weight section coordinate closest to p
	static double meas_section(const Point& p, const Point& L1, const Point& L2) noexcept;
	static double meas_section(const Point& p, const Point& L1, const Point& L2, double& ksi) noexcept;
	//measure from p to [L1, L2] line
	static double meas_line(const Point& p, const Point& L1, const Point& L2) noexcept;
	//normalized line equation
	static std::array<double, 3> line_eq(const Point& p1, const Point& p2) noexcept;

	//finds a point between two given ones with a certain weight w
	static Point Weigh(const Point& p1, const Point& p2, double w) noexcept{
		return Point(p1.x*(1-w)+p2.x*w, p1.y*(1-w)+p2.y*w);
	}


	//finds point which coordinates are the highest among all points in container
	template<class FirstIter, class LastIter>
	static Point GetTop(FirstIter start, LastIter end){
		return Point(
			std::max_element(start, end, 
				[](const Point& p1, const Point& p2){ return p1.x<p2.x; })->x,
			std::max_element(start, end, 
				[](const Point& p1, const Point& p2){ return p1.y<p2.y; })->y
		);
	}
	template<class FirstIter, class LastIter>
	static Point GetBot(FirstIter start, LastIter end){
		return Point(
			std::min_element(start, end, 
				[](const Point& p1, const Point& p2){ return p1.x<p2.x; })->x,
			std::min_element(start, end, 
				[](const Point& p1, const Point& p2){ return p1.y<p2.y; })->y
		);
	}

	template<class FirstIter, class LastIter>
	static Point CPoint(FirstIter start, LastIter end){
		return std::accumulate(start, end, Point(0,0))/std::distance(start,end);
	}

};

//point operators
inline Point operator+(const Point& left, const Point& right) noexcept{
	return Point(left.x+right.x, left.y+right.y);
}
inline Point operator-(const Point& left, const Point& right) noexcept{
	return Point(left.x-right.x, left.y-right.y);
}
inline Point& operator+=(Point& left, const Point& right) noexcept{
	left.x+=right.x; left.y+=right.y;
	return left;
}
inline Point& operator-=(Point& left, const Point& right) noexcept{
	left.x-=right.x; left.y-=right.y;
	return left;
}
inline Point& operator/=(Point& left, double d) noexcept{
	left.x/=d; left.y/=d;
	return left;
}
inline Point& operator*=(Point& left, double d) noexcept{
	left.x*=d; left.y*=d;
	return left;
}
inline Point operator/(const Point& p, double d) noexcept{ auto x=Point(p); return std::move(x/=d); }
inline Point operator*(const Point& p, double d) noexcept{ auto x=Point(p); return std::move(x*=d); }
inline bool operator==(const Point& p1, const Point& p2) noexcept{ return (ISEQ(p1.x, p2.x) && ISEQ(p1.y, p2.y)); }
inline bool operator!=(const Point& p1, const Point& p2) noexcept{ return (!ISEQ(p1.x, p2.x) || !ISEQ(p1.y, p2.y)); }
inline bool operator<(const Point& p1, const Point& p2) noexcept{
	return (ISEQ(p1.x, p2.x)) ? ISLOWER(p1.y, p2.y) : ISLOWER(p1.x, p2.x);
}

inline std::ostream& operator<<(std::ostream& os, const Point& p){
	os<<p.x<<" "<<p.y;
	return os;
}

bool isOnSection(const Point& p, const Point& start, const Point& end, double& ksi, double eps=geps) noexcept;

//Finds a cross point between two sections: (p1S,p1E) and (p2S, p2E).
//ksieta -- ouput cross weights: [0] -- weight in first section, [1] -- weight in second section
//Returns true if 0<=ksieta[0,1]<=1.
//if sections are parallel: ksieta[0,1]=gbig, returns false
bool SectCross(const Point& p1S, const Point& p1E, const Point& p2S, const Point& p2E, double* ksieta) noexcept;
//same but with internal renorming. Gives better results aber it is slower.
bool SectCrossWRenorm(const Point& p1S, const Point& p1E, const Point& p2S, const Point& p2E, double* ksieta) noexcept;

//=>
//  0 if p lies to the left of [L1->L2] line
//  1 if p lies on [L1->L2] line
//  2 if p lies to the right of [L1->L2] line
// -1 if (L1 == L2) 
int LinePointWhereIs(const Point& p, const Point& L1, const Point& L2) noexcept;

//scaling
struct ScaleBase{
	Point p0;
	double L;
	explicit ScaleBase(double x0=0, double y0=0, double _L=1): p0(x0, y0), L(_L){}
	//scale and unscale procedures
	void scale(Point& p) const noexcept {p-=p0; p/=L; }
	void unscale(Point& p) const noexcept {p*=L; p+=p0; }
	//scale points in container with Points
	template<class FirstIter, class LastIter>
	void scale(FirstIter start, LastIter end) const noexcept {
		for (auto it=start; it!=end; ++it){ scale(*it); }
	}
	template<class FirstIter, class LastIter>
	void unscale(FirstIter start, LastIter end) const noexcept{
		for (auto it=start; it!=end; ++it){ unscale(*it); }
	}
	//forces all points in container be in 1x1x1 square preserving aspect ratio
	template <class Container>
	static ScaleBase doscale(Container& cont) noexcept{
		return doscale(cont.begin(), cont.end());
	}
	template<class FirstIter, class LastIter>
	static ScaleBase doscale(FirstIter start, LastIter end) noexcept{
		if (start >= end) return ScaleBase();
		auto p0 = Point::GetBot(start, end);
		auto pdif = Point::GetTop(start, end) - p0;
		double L = std::max(pdif.x, pdif.y);
		ScaleBase ret(p0.x, p0.y, L);
		ret.scale(start, end);
		return ret;
	}
	template<class Iter>
	static ScaleBase p_doscale(Iter start, Iter end, double a=1.0) noexcept{
		if (start >= end) return ScaleBase();
		vector<Point> p; p.reserve(end-start);
		std::transform(start, end, std::back_inserter(p), 
				[](const typename Iter::value_type& x){ return *x; });
		auto p0 = Point::GetBot(p.begin(), p.end());
		auto pdif = Point::GetTop(p.begin(), p.end()) - p0;
		double L = std::max(pdif.x, pdif.y) / a;
		ScaleBase ret(p0.x, p0.y, L);
		ret.p_scale(start, end);
		return ret;
	}
	//scale points in container with Points pointers
	template<class FirstIter, class LastIter>
	void p_scale(FirstIter start, LastIter end) const noexcept {
		for (auto it=start; it!=end; ++it){ scale(**it); }
	}
	template<class FirstIter, class LastIter>
	void p_unscale(FirstIter start, LastIter end) const noexcept{
		for (auto it=start; it!=end; ++it){ unscale(**it); }
	}
};

//find a node by coordinate
class NodeFinder{
	const int Nx, Ny;
	const double eps;
	const double x0, y0, x1, y1, hx, hy;
	vector<int> get_index(const Point* p) const;
	int to_glob(int i, int j) const { return j*Nx + i; }
	std::vector<std::vector<const Point*>> data;
	bool is_equal_point(const Point* p1, const Point* p2) const;
public:
	NodeFinder(Point p0, double Lx, double Ly, int Nx=30, int Ny=30, double eps=geps);
	NodeFinder(const std::pair<Point, Point>& rect, int Nx=30, int Ny=30, double eps=geps);
	//adds point to data list if necessary
	//returns pointer to previously added point if p was already added,
	//p if point was added to point data
	//and 0 if p lies outside defined rectangle
	const Point* add(const Point* p);
	//returns pointer to previously added point if p was already added
	//and 0 otherwise
	const Point* find(const Point* p) const;
};

//Angles
//to angle in [0, 2pi]
inline double ToAngle(double angle, double eps=0.0){
	if (fabs(angle)<eps || fabs(angle-2*M_PI)<eps) return 0.0;
	if (angle<0) return ToAngle(angle+2*M_PI);
	else if (angle>=2*M_PI) return ToAngle(angle-2*M_PI, eps);
	else return angle;
}
inline double DegToAngle(double deg){
	return ToAngle(deg/180*M_PI);
}
inline double AngleAdd(double angle, double add, double eps=0.0){ return ToAngle(angle+add,eps); }
inline double Angle(const Point& p1, const Point& p2, const Point& p3){
	return ToAngle(atan2(p1.y-p2.y, p1.x-p2.x) - atan2(p3.y-p2.y, p3.x-p2.x));
}

//Vector
typedef Point Vect;
inline double vecDot(const Vect& u, const Vect& v){ return u.x*v.x+u.y*v.y; }
inline double vecLen(const Vect& a){ return sqrt(vecDot(a,a)); }
inline void vecSetLen(Vect& a, double Len){ a/=(vecLen(a)/Len); }
inline void vecNormalize(Vect& a){ a/=vecLen(a); }
inline double vecCrossZ(const Vect& a, const Vect& b){ return a.x*b.y-a.y*b.x; }
inline Vect vecRotate(const Vect& u, double angle){
	double cs = cos(angle), sn = sin(angle);
	return Point(u.x*cs - u.y*sn, u.x*sn+u.y*cs);
}

//triangle procedures
inline double triarea(const Point& p1, const Point& p2, const Point& p3){
	return 0.5*vecCrossZ(p1-p3, p2-p3);
}

//add refinement points within [0, Len] section. 0 and Len are not included.
vector<double> RefineSection(double a, double b, double Len, double Den);

class BoundingBox{
protected:
	void init();
public:
	double xmin, xmax, ymin, ymax;

	BoundingBox():xmin(0), xmax(1), ymin(0), ymax(1){}
	BoundingBox(double x0, double y0, double x1, double y1):xmin(x0), xmax(x1), ymin(y0), ymax(y1){}
	BoundingBox(const vector<BoundingBox>&, double e=0.0);
	BoundingBox(const Point& p1, const Point& p2, double e=0.0);
	BoundingBox(const Point& p1, double e=0.0);

	//Container<Point> -> BoundingBox
	template<class Iter> static typename std::enable_if<
		std::is_base_of<Point, typename Iter::value_type>::value,
		BoundingBox
	>::type Build(Iter first, Iter last){
		BoundingBox ret(first->x, first->y, first->x, first->y);
		while (first!=last) {ret.WidenWithPoint(&*first); ++first;}
		return ret;
	}
	//Container<Point*> -> BoundingBox
	template<class Iter> static typename std::enable_if<
		!std::is_base_of<Point, typename Iter::value_type>::value,
		BoundingBox
	>::type Build(Iter first, Iter last){
		BoundingBox ret((*first)->x, (*first)->y, (*first)->x, (*first)->y);
		while (first!=last) {ret.WidenWithPoint(**first); ++first;}
		return ret;
	}

	//Enlarge if point lies outside box
	void WidenWithPoint(const Point& p);
	//widen to all directions at certain distance
	void widen(double e);

	Point Center() const;

	double area() const;
	double lenx() const { return xmax-xmin; }
	double leny() const { return ymax-ymin; }
	double maxlen() const { return std::max(lenx(), leny()); }
	double lendiag() const { return sqrt(lenx()*lenx() + leny()*leny()); }
	Point BottomLeft() const { return Point(xmin, ymin); }
	Point TopRight() const { return Point(xmax, ymax); }
	vector<Point> FourPoints() const;
	ScaleBase to_scale() const { return ScaleBase(xmin, ymin, std::max(lenx(), leny())); }

	//-> INSIDE, OUTSIDE, BOUND
	int whereis(const Point& p) const;
	//does this have any intersections or tangent segments with another segment
	bool has_common_points(const BoundingBox& bb) const;
	//does this contain any part of [p1, p2] segment
	bool contains(const Point& p1, const Point& p2) const;

	//Filter points from container of Point*
	template<class Container>
	Container Filter(const Container& data,
			bool inside, bool bound, bool outside){
		Container out;
		for (auto& p: data){
			int pos = whereis(*p);
			if (inside && pos == INSIDE) out.push_back(p);
			else if (bound && pos == BOUND) out.push_back(p);
			else if (outside && pos == OUTSIDE) out.push_back(p);
		}
		return out;
	}
};

struct BoundingBoxFinder{
	BoundingBoxFinder(const BoundingBox& area, double L);

	void addentry(const BoundingBox& e);
	vector<int> suspects(const BoundingBox& bb) const;
	vector<int> suspects(const Point& bb) const;
private:
	double x0, y0, hx, hy;
	int mx, my;
	int totaldata;
	vector<vector<int>> data;
	vector<int> allsuspects(int x0, int x1, int y0, int y1) const;
	void addsuspects(int ix, int iy, vector<int>& ret) const;
	int get_xstart(double x) const;
	int get_ystart(double y) const;
	int get_xend(double x) const;
	int get_yend(double y) const;
};

#endif
