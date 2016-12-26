#include "bgeom3d.h"
#include "hmcompute.hpp"

double Point3::sdist_plane(const Point3& p, const Point3& plane1, const Point3& plane2, const Point3& plane3){
	Vect3 n = right_normal(plane1, plane2, plane3);
	return vecDot(n, p-plane1);
}

namespace{

int pos1D(double x, double a, double b){
	if (ISLOWER(x, a)) return -1;
	else if (ISGREATER(x, b)) return 1;
	else return 0;
			
};
}
//returns INSIDE, OUTSIDE, BOUND
int BoundingBox3D::position(const BoundingBox3D& bb) const{
	int x1 = pos1D(bb.xmin, xmin, xmax);
	int x2 = pos1D(bb.xmax, xmin, xmax);
	if (x1 == x2 && x1 != 0) return OUTSIDE;
	int y1 = pos1D(bb.ymin, ymin, ymax);
	int y2 = pos1D(bb.ymax, ymin, ymax);
	if (y1 == y2 && y1 != 0) return OUTSIDE;
	int z1 = pos1D(bb.zmin, zmin, zmax);
	int z2 = pos1D(bb.zmax, zmin, zmax);
	if (z1 == z2 && z1 != 0) return OUTSIDE;

	if (!x1 && !x2 && !y1 && !y2 && !z1 && !z2) return INSIDE;
	else if ((x1 == x2 && x1 != 0) || (y1 == y2 && y1 != 0) || (z1 == z2 && z1 != 0)) return OUTSIDE;
	else return BOUND;
}
int BoundingBox3D::position(const Point3& v) const{
	if (ISLOWER(v.x, xmin) || ISGREATER(v.x, xmax) ||
		ISLOWER(v.y, ymin) || ISGREATER(v.y, ymax) ||
		ISLOWER(v.z, zmin) || ISGREATER(v.z, zmax)) return OUTSIDE;
	else if (ISGREATER(v.x, xmin) && ISLOWER(v.x, xmax) &&
		ISGREATER(v.y, ymin) && ISLOWER(v.y, ymax) &&
		ISGREATER(v.z, zmin) && ISLOWER(v.z, zmax)) return INSIDE;
	return BOUND;
}
//returns:
//0 - bb equals this
//1 - bb is inside this (with possible touched faces)
//2 - this is inside bb (with possible touched faces)
//3 - bb and this doesn't intersect (with possible touched faces)
//4 - faces of bb and this cross each other
int BoundingBox3D::relation(const BoundingBox3D& bb) const{
	auto pos1Dex = [](double x, double a, double b)->int{
		if (ISLOWER(x, a)) return -2;
		if (ISGREATER(x, b)) return 2;
		if (ISEQ(x, a)) return -1;
		if (ISEQ(x, b)) return 1;
		return 0;
	};
	int x1 = pos1Dex(xmin, bb.xmin, bb.xmax);
	if (x1 >= 1) return 3;
	int x2 = pos1Dex(xmax, bb.xmin, bb.xmax);
	if (x2 <= -1) return 3;
	int y1 = pos1Dex(ymin, bb.ymin, bb.ymax);
	if (y1 >= 1) return 3;
	int y2 = pos1Dex(ymax, bb.ymin, bb.ymax);
	if (y2 <= -1) return 3;
	int z1 = pos1Dex(zmin, bb.zmin, bb.zmax);
	if (z1 >= 1) return 3;
	int z2 = pos1Dex(zmax, bb.zmin, bb.zmax);
	if (z2 <= -1) return 3;

	if (x1 == -1  && x2 == 1 && y1 == -1 && y2 == 1 && z1 == -1 && z2 == 1) return 0;
	if (x1 <= -1 && x2 >= 1 && y1 <= -1 && y2 >= 1 && z1 <= -1 && z2 >= 1) return 1;
	if (x1 >= -1 && x2 <=1 && y1 >= -1 && y2 <= 1 && z1 >= -1 && z2 <= 1) return 2;
	return 4;
}

double tetrahedron_volume_0(const Point3& p1, const Point3& p2, const Point3& p3){
	double A[] = {p1.x, p1.y, p1.z,
		      p2.x, p2.y, p2.z,
		      p3.x, p3.y, p3.z};
	return HMMath::Compute::determinant_3x3(A)/6.;
}
