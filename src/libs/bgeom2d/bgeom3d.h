#ifndef HYBMESH_BGEOM3D_H
#define HYBMESH_BGEOM3D_H
#include "hmproject.h"

struct Point3{
	double x,y,z;
	Point3(double x=0, double y=0, double z=0):x(x), y(y), z(z){}
	Point3& operator=(const Point3& p){
		if (this!=&p){
			this->x=p.x; this->y=p.y;
			this->z=p.z;
		}
		return *this;
	}
	void set(double a, double b, double c){x=a; y=b; z=c;}
	void set(const Point3& p2){x=p2.x; y=p2.y; z=p2.z;}

	static double meas(const Point3& p1, const Point3& p2){
		return (p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y)+(p1.z-p2.z)*(p1.z-p2.z);
	}
	static double dist(const Point3& p1, const Point3& p2){
		return sqrt(meas(p1,p2));
	}
};
//point operators
inline Point3 operator+(const Point3& left, const Point3& right) {
	return Point3(left.x+right.x, left.y+right.y, left.z+right.z);
}
inline Point3 operator-(const Point3& left, const Point3& right) {
	return Point3(left.x-right.x, left.y-right.y, left.z-right.z);
}
inline Point3& operator+=(Point3& left, const Point3& right) {
	left.x+=right.x; left.y+=right.y; left.z+=right.z;
	return left;
}
inline Point3& operator-=(Point3& left, const Point3& right) {
	left.x-=right.x; left.y-=right.y; left.z-=right.z;
	return left;
}
inline Point3& operator/=(Point3& left, double d) {
	left.x/=d; left.y/=d; left.z/=d;
	return left;
}
inline Point3& operator*=(Point3& left, double d) {
	left.x*=d; left.y*=d; left.z*=d;
	return left;
}
inline Point3 operator/(const Point3& p, double d) { auto x=Point3(p); return std::move(x/=d); }
inline Point3 operator*(const Point3& p, double d) { auto x=Point3(p); return std::move(x*=d); }
inline bool operator==(const Point3& p1, const Point3& p2) { return (ISEQ(p1.x, p2.x) && ISEQ(p1.y, p2.y) && ISEQ(p1.z, p2.z)); }
inline bool operator!=(const Point3& p1, const Point3& p2) { return (!ISEQ(p1.x, p2.x) || !ISEQ(p1.y, p2.y) || !ISEQ(p1.z, p2.z)); }
inline bool operator<(const Point3& p1, const Point3& p2) {
	if (ISEQ(p1.x, p2.x)){
		if (ISEQ(p1.y, p2.y)){
			return ISLOWER(p1.z, p2.z);
		} else return ISLOWER(p1.y, p2.y);
	} else return ISLOWER(p1.x, p2.x);
}
inline std::ostream& operator<<(std::ostream& os, const Point3& p){
	os<<p.x<<" "<<p.y<<" "<<p.z<<std::endl;
	return os;
}

struct BoundingBox3D{
	double xmin, xmax, ymin, ymax, zmin, zmax;

	BoundingBox3D(): xmin(0), xmax(0), ymin(0), ymax(0), zmin(0), zmax(0){};
	template<class PContainer>
	BoundingBox3D(const PContainer& p){
		if (p.size() == 0){
			xmin = xmax = ymin = ymax = zmin = zmax = 0;
		} else {
			auto it = p.begin();
			xmin = xmax = (*it)->x;
			ymin = ymax = (*it)->y;
			zmin = zmax = (*it)->z;
			++it;
			while (it != p.end()){
				if ((*it)->x > xmax) xmax = (*it)->x;
				if ((*it)->x < xmin) xmin = (*it)->x;
				if ((*it)->y > ymax) ymax = (*it)->y;
				if ((*it)->y < ymin) ymin = (*it)->y;
				if ((*it)->z > zmax) zmax = (*it)->z;
				if ((*it)->z < zmin) zmin = (*it)->z;
				++it;
			}
		}
	}

	//returns INSIDE if object is fully inside this,
	//        OUTSIDE if there is no intersection,
	//        BOUND otherwise
	int position(const BoundingBox3D& bb) const;
	int position(const Point3& v) const;

	//returns:
	//0 - bb equals this
	//1 - bb is inside this (with possible touched faces)
	//2 - this is inside bb (with possible touched faces)
	//3 - bb and this doesn't intersect (with possible touched faces)
	//4 - faces of bb and this cross each other
	int relation(const BoundingBox3D& bb) const;

	BoundingBox3D widen(double eps){
		xmin-=eps; xmax+=eps;
		ymin-=eps; ymax+=eps;
		zmin-=eps; zmax+=eps;
	}

	double xlen() const { return xmax - xmin; }
	double ylen() const { return ymax - ymin; }
	double zlen() const { return zmax - zmin; }
	double maxlen() const { return std::max(xlen(), std::max(ylen(), zlen())); }
	Point3 center() const { return Point3((xmax+xmin)/2, (ymax+ymin)/2, (zmax+zmin)/2); }

};

struct ScaleBase3{
	Point3 p0;
	double L;
	ScaleBase3(double x0=0, double y0=0, double z0=0, double L=1): p0(x0, y0, z0), L(L){}
	void scale(Point3& p) const  {p-=p0; p/=L; }
	void unscale(Point3& p) const  {p*=L; p+=p0; }
	//scale and unscale procedures
	//inscribes all points into [0, a]x[0, a]x[0, a] cube
	template<class PContainer>
	static ScaleBase3 doscale(PContainer& p, double a=1.0) {
		if (p.size() == 0) return ScaleBase3();
		BoundingBox3D bb(p);
		ScaleBase3 ret(bb.xmin, bb.ymin, bb.zmin, bb.maxlen()/a);
		ret.scale(p.begin(), p.end());
		return ret;
	}

	//scale points in container with Points pointers
	template<class FirstIter, class LastIter>
	void scale(FirstIter begin, LastIter end) const  {
		for (auto it=begin; it!=end; ++it){ scale(**it); }
	}
	template<class FirstIter, class LastIter>
	void unscale(FirstIter begin, LastIter end) const  {
		for (auto it=begin; it!=end; ++it){ unscale(**it); }
	}
};
//signed tetrahedron volume which first point is (0, 0, 0)
double tetrahedron_volume_0(const Point3& p1, const Point3& p2, const Point3& p3);
//signed tetrahedron volume by four points
inline double tetrahedron_volume(const Point3& p0, const Point3& p1, const Point3& p2, const Point3& p3){
	return tetrahedron_volume_0(p1-p0, p2-p0, p3-p0);
}


typedef Point3 Vect3;
inline double vecDot(const Vect3& a, const Vect3& b){ return (a.x*b.x + a.y*b.y + a.z*b.z); }
inline double vecLen(const Vect3& a){ return sqrt(vecDot(a, a)); }
inline Vect3 vecCross(const Vect3& a, const Vect3& b){
	return Vect3(a.y*b.z-a.z*b.y,
	            -a.x*b.z+a.z*b.x,
		     a.x*b.y-a.y*b.x);
}
inline void vecNormalize(Vect3& a){
	double L = vecLen(a);
	a/=L;
}

#endif
