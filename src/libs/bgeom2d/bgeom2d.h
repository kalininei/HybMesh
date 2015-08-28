#ifndef HYBMESH_BGEOM_H
#define HYBMESH_BGEOM_H

#include "hmproject.h"
#include <algorithm>
#include "addalgo.hpp"

const double geps = 1e-8;
const double geps2 = 1e-14;
const double gbig = 1.0/geps;
inline bool ISZERO(double v){ return fabs(v)<geps; }
inline bool ISEQ(double a, double b){ return ISZERO(a-b); }
inline bool ISEQLOWER(double x, double y){ return x<y+geps; }
inline bool ISEQGREATER(double x, double y){ return x>y-geps; }
inline bool ISLOWER(double x, double y){ return x<y-geps; }
inline bool ISGREATER(double x, double y){ return x>y+geps; }
inline int SIGN(double x){ if (ISZERO(x)) return 0; else return (x<0) ? -1 : 1; }
inline double sqr(double x){ return x*x; }

//positioning constants
const int BOUND = 0;
const int INSIDE = 1;
const int OUTSIDE = -1;

//Point
struct Point{
	double x,y;
	//constructors
	explicit Point(double _x=0, double _y=0):x(_x), y(_y){}
	Point(const Point& p2):x(p2.x), y(p2.y){}
	Point& operator=(const Point& p){
		if (this!=&p){
			this->x=p.x; this->y=p.y;
		}
		return *this;
	}
	void set(double a=0, double b=0) noexcept {x=a; y=b;}
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
	static double meas_section(const Point& p, const Point& L1, const Point& L2) noexcept;

	//finds a point between two given ones with a certain weight w
	static Point Weigh(const Point& p1, const Point& p2, double w){
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
inline Point operator+(const Point& left, const Point& right){
	return Point(left.x+right.x, left.y+right.y);
}
inline Point operator-(const Point& left, const Point& right){
	return Point(left.x-right.x, left.y-right.y);
}
inline Point& operator+=(Point& left, const Point& right){
	left.x+=right.x; left.y+=right.y;
	return left;
}
inline Point& operator-=(Point& left, const Point& right){
	left.x-=right.x; left.y-=right.y;
	return left;
}
inline Point& operator/=(Point& left, double d){
	left.x/=d; left.y/=d;
	return left;
}
inline Point& operator*=(Point& left, double d){
	left.x*=d; left.y*=d;
	return left;
}
inline Point operator/(const Point& p, double d){ auto x=Point(p); return std::move(x/=d); }
inline Point operator*(const Point& p, double d){ auto x=Point(p); return std::move(x*=d); }
inline bool operator==(const Point& p1, const Point& p2){ return (ISEQ(p1.x, p2.x) && ISEQ(p1.y, p2.y)); }
inline bool operator<(const Point& p1, const Point& p2){
	return (ISEQ(p1.x, p2.x)) ? ISLOWER(p1.y, p2.y) : ISLOWER(p1.x, p2.x);
}

inline std::ostream& operator<<(std::ostream& os, const Point& p){
	os<<p.x<<" "<<p.y;
	return os;
}

bool isOnSection(const Point& p, const Point& start, const Point& end, double& ksi, double eps=geps);

//Finds a cross point between two sections: (p1S,p1E) and (p2S, p2E).
//ksieta -- ouput cross weights: [0] -- weight in first section, [1] -- weight in second section
//Returns true if 0<=ksieta[0,1]<=1.
//if sections are parallel: ksieta[0,1]=gbig, returns false
bool SectCross(const Point& p1S, const Point& p1E, const Point& p2S, const Point& p2E, double* ksieta) noexcept;

//scaling
struct ScaleBase{
	const Point p0;
	const double L;
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
	template<class FirstIter, class LastIter>
	static ScaleBase doscale(FirstIter start, LastIter end) noexcept{
		auto p0 = Point::GetBot(start, end);
		auto pdif = Point::GetTop(start, end) - p0;
		double L = std::max(pdif.x, pdif.y);
		ScaleBase ret(p0.x, p0.y, L);
		ret.scale(start, end);
		return ret;
	}
	template<class Iter>
	static ScaleBase p_doscale(Iter start, Iter end) noexcept{
		vector<Point> p; p.reserve(end-start);
		std::transform(start, end, std::back_inserter(p), 
				[](const typename Iter::value_type& x){ return *x; });
		auto p0 = Point::GetBot(p.begin(), p.end());
		auto pdif = Point::GetTop(p.begin(), p.end()) - p0;
		double L = std::max(pdif.x, pdif.y);
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
inline double ToAngle(double angle, double eps=0.0){
	if (fabs(angle)<eps || fabs(angle-2*M_PI)<eps) return 0.0;
	if (angle<0) return ToAngle(angle+2*M_PI);
	else if (angle>=2*M_PI) return ToAngle(angle-2*M_PI, eps);
	else return angle;
}
inline double AngleAdd(double angle, double add, double eps=0.0){ return ToAngle(angle+add,eps); }

//Vector
typedef Point Vect;
inline double vecDot(const Vect& u, const Vect& v){ return u.x*v.x+u.y*v.y; }
inline double vecLen(const Vect& a){ return sqrt(vecDot(a,a)); }
inline void vecSetLen(Vect& a, double Len){ a/=(vecLen(a)/Len); }
inline void vecNormalize(Vect& a){ a/=vecLen(a); }
inline double vecCrossZ(const Vect& a, const Vect& b){ return a.x*b.y-a.y*b.x; }

//triangle procedures
inline double triarea(const Point& p1, const Point& p2, const Point& p3){
	return 0.5*vecCrossZ(p1-p3, p2-p3);
}

//add refinement points within [0, Len] section. 0 and Len are not included.
vector<double> RefineSection(double a, double b, double Len, double Den);

class BoundingBox{
protected:
	void init();
	void add_point(const Point* p);
	void widen(double e);
public:
	double xmin, xmax, ymin, ymax;

	BoundingBox():xmin(0), xmax(1), ymin(0), ymax(1){}
	BoundingBox(double x0, double y0, double x1, double y1):xmin(x0), xmax(x1), ymin(y0), ymax(y1){}
	BoundingBox(const vector<BoundingBox>&, double e=0.0);
	BoundingBox(const Point& p1, const Point& p2, double e=0.0);

	void WidenWithPoint(const Point& p);

	double area() const;
	double lenx() const { return xmax-xmin; }
	double leny() const { return ymax-ymin; }
	Point BottomLeft() const { return Point(xmin, ymin); }
	Point TopRight() const { return Point(xmax, ymax); }

	//-> INSIDE, OUTSIDE, BOUND
	int whereis(const Point& p) const;
	//does this have any intersections or tangent segments with another segment
	bool has_common_points(const BoundingBox& bb) const;
	//does this contain any part of [p1, p2] segment
	bool contains(const Point& p1, const Point& p2) const;
};

#endif
