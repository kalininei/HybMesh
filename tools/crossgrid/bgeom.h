#ifndef CROSSGRID_BGEOM_H
#define CROSSGRID_BGEOM_H

#include <math.h>
#include <algorithm>
#include <iostream>
#include <set>
#include <memory>
#include "addalgo.hpp"

using std::vector;
template<class T> using shp_vector = std::vector<std::shared_ptr<T>>;

const double geps = 1e-8;
const double geps2 = geps*geps;
const double gbig = 1.0/geps;
inline bool ISZERO(double v){ return fabs(v)<geps; }
inline bool ISEQ(double a, double b){ return ISZERO(a-b); }
inline bool ISEQLOWER(double x, double y){ return x<y+geps; }
inline bool ISEQGREATER(double x, double y){ return x>y-geps; }
inline bool ISLOWER(double x, double y){ return x<y-geps; }
inline bool ISGREATER(double x, double y){ return x>y+geps; }

inline double sqr(double x){ return x*x; }

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
			double x=(Dim>0) ? *it++ : 0.0; 
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
	//scale points in container with Points pointers
	template<class FirstIter, class LastIter>
	void p_scale(FirstIter start, LastIter end) const noexcept {
		for (auto it=start; it!=end; ++it){ scale(**it); }
	}
	template<class FirstIter, class LastIter>
	void p_unscale(FirstIter start, LastIter end) const noexcept{
		for (auto it=start; it!=end; ++it){ unscale(**it); }
	}
	//forces all points in container be in 1x1x1 square preserving aspect ratio
	template<class Iter>
	static ScaleBase p_doscale(Iter start, Iter end) noexcept{
		vector<Point> p; p.reserve(end-start);
		std::transform(start, end, std::back_inserter(p), 
				[](const typename Iter::value_type& x){ return *x; });
		return doscale(p.begin(), p.end());
	}
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


#endif
