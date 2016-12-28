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
	
	//finds a point between two given ones with a certain weight w
	static Point3 Weigh(const Point3& p1, const Point3& p2, double w){
		return Point3((1-w)*p1.x + w*p2.x,
			(1-w)*p1.y + w*p2.y,
			(1-w)*p1.z + w*p2.z);
	}

	static double meas(const Point3& p1, const Point3& p2){
		return (p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y)+(p1.z-p2.z)*(p1.z-p2.z);
	}
	static double dist(const Point3& p1, const Point3& p2){
		return sqrt(meas(p1,p2));
	}
	//signed distance to plane. Gives positive if p is to the right of the plane
	static double sdist_plane(const Point3& p, const Point3& plane1, const Point3& plane2, const Point3& plane3);
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
inline Point3 operator/(const Point3& p, double d) { auto x=Point3(p); return x/=d; }
inline Point3 operator*(const Point3& p, double d) { auto x=Point3(p); return x*=d; }
inline Point3 operator*(double d, const Point3& p) { auto x=Point3(p); return x*=d; }
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

//local coordinates of intersection between triangle (tri1, tri2, tri3) and segment (pstart, pend)
//calculates: xke[0] - local coordinate on (pstart-pend)
//            xke[1] - local plane coordinate by (tri2-tri1) basis
//            xke[2] - local plane coordinate by (tri3-tri1) basis
//returns true if all local coordinates are in [0, 1] and (xke[1] + xke[2] in [0, 1])
//    using equal or greater comparison.
//if segment and triangle are parallel then returns false and xke[:]=gbig
bool segment_triangle_cross3d(const Point3& pstart, const Point3& pend,
		const Point3& tri1, const Point3& tri2, const Point3& tri3,
		double* xke);

struct BoundingBox3D{
	double xmin, xmax, ymin, ymax, zmin, zmax;

	BoundingBox3D(): xmin(0), xmax(0), ymin(0), ymax(0), zmin(0), zmax(0){};
	BoundingBox3D(const Point3& p1, const Point3& p2){
		xmin = p1.x; xmax = p2.x;
		ymin = p1.y; ymax = p2.y;
		zmin = p1.z; zmax = p2.z;
		if (xmin > xmax) std::swap(xmin, xmax);
		if (ymin > ymax) std::swap(ymin, ymax);
		if (zmin > zmax) std::swap(zmin, zmax);
	}
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

	void widen(double eps){
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

struct BoundingBox3DFinder{
	//L - step size
	BoundingBox3DFinder(const BoundingBox3D& area, double L);

	void addentry(const BoundingBox3D& e);
	vector<int> suspects(const BoundingBox3D& bb) const;
	vector<int> suspects(const Point3& bb) const;
private:
	double x0, y0, z0, hx, hy, hz;
	int mx, my, mz;
	int totaldata;
	vector<vector<int>> data;
	vector<int> allsuspects(int x0, int x1, int y0, int y1, int z0, int z1) const;
	void addsuspects(int ix, int iy, int iz, vector<int>& ret) const;
	int get_xstart(double x) const;
	int get_ystart(double y) const;
	int get_zstart(double y) const;
	int get_xend(double x) const;
	int get_yend(double y) const;
	int get_zend(double y) const;
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
inline Vect3 right_normal_0(const Point3& p2, const Point3& p3){
	Vect3 ret = vecCross(p2, p3);
	vecNormalize(ret);
	return ret;
}
inline Vect3 left_normal_0(const Point3& p2, const Point3& p3){
	Vect3 ret = vecCross(p3, p2);
	vecNormalize(ret);
	return ret;
}
inline Vect3 right_normal(const Point3& p1, const Point3& p2, const Point3& p3){
	return right_normal_0(p2-p1, p3-p1);
}
inline Vect3 left_normal(const Point3& p1, const Point3& p2, const Point3& p3){
	return left_normal_0(p2-p1, p3-p1);
}

#endif
