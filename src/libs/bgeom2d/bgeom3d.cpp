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

bool segment_triangle_cross3d(const Point3& pstart, const Point3& pend,
		const Point3& tri1, const Point3& tri2, const Point3& tri3,
		double* xke){
	double A[9]={
		pstart.x-pend.x, tri2.x-tri1.x, tri3.x-tri1.x,
		pstart.y-pend.y, tri2.y-tri1.y, tri3.y-tri1.y,
		pstart.z-pend.z, tri2.z-tri1.z, tri3.z-tri1.z};
	double rhs[3]={
		pstart.x - tri1.x,
		pstart.y - tri1.y,
		pstart.z - tri1.z};
	if (!HMMath::Compute::mat_solve_3x3(A, rhs, xke)){
		xke[0] = gbig; xke[1] = gbig; xke[2] = gbig;
		return false;
	} else {
		return ISIN_EE(xke[0], 0, 1) && ISIN_EE(xke[1], 0, 1)
			&& ISIN_EE(xke[2], 0, 1-xke[1]);
	}
}

BoundingBox3DFinder::BoundingBox3DFinder(const BoundingBox3D& area, double L){
	totaldata = 0;
	x0 = area.xmin;
	y0 = area.ymin;
	z0 = area.zmin;
	mx = std::ceil(area.xlen()/L)+1;
	my = std::ceil(area.ylen()/L)+1;
	mz = std::ceil(area.zlen()/L)+1;
	hx = area.xlen()/(mx-1);
	hy = area.ylen()/(my-1);
	hz = area.zlen()/(mz-1);
	data.resize((mx+1)*(my+1)*(mz+1));
}
void BoundingBox3DFinder::addentry(const BoundingBox3D& e){
	int x0 = get_xstart(e.xmin);
	int x1 = get_xend(e.xmax);
	int y0 = get_ystart(e.ymin);
	int y1 = get_yend(e.ymax);
	int z0 = get_zstart(e.zmin);
	int z1 = get_zend(e.zmax);
	for (int i=x0; i<=x1; ++i)
	for (int j=y0; j<=y1; ++j)
	for (int k=z0; k<=z1; ++k){
		int gi = k*(mx+1)*(my+1)+j*(mx+1)+i;
		data[gi].push_back(totaldata);
	}
	++totaldata;
}
int BoundingBox3DFinder::get_xstart(double x) const{
	int ret = (x - x0) / hx;
	if (ISEQ(ret*hx, x-x0)) --ret;
	if (ret<0) ret = 0;
	return ret;
}
int BoundingBox3DFinder::get_ystart(double y) const{
	int ret = (y - y0) / hy;
	if (ISEQ(ret*hy, y-y0)) --ret;
	if (ret<0) ret = 0;
	return ret;
}
int BoundingBox3DFinder::get_zstart(double z) const{
	int ret = (z - z0) / hz;
	if (ISEQ(ret*hz, z-z0)) --ret;
	if (ret<0) ret = 0;
	return ret;
}
int BoundingBox3DFinder::get_xend(double x) const{
	int ret = (x - x0) / hx;
	if (ISEQ((ret+1)*hx, x-x0)) ++ret;
	if (ret>mx) ret = mx;
	return ret;
}
int BoundingBox3DFinder::get_yend(double y) const{
	int ret = (y - y0) / hy;
	if (ISEQ((ret+1)*hy, y-y0)) ++ret;
	if (ret>my) ret = my;
	return ret;
}
int BoundingBox3DFinder::get_zend(double z) const{
	int ret = (z - z0) / hz;
	if (ISEQ((ret+1)*hz, z-z0)) ++ret;
	if (ret>mz) ret = mz;
	return ret;
}

vector<int> BoundingBox3DFinder::suspects(const BoundingBox3D& bb) const{
	return allsuspects(
		get_xstart(bb.xmin), get_xend(bb.xmax),
		get_ystart(bb.ymin), get_yend(bb.ymax),
		get_zstart(bb.zmin), get_zend(bb.zmax));
}
vector<int> BoundingBox3DFinder::suspects(const Point3& bb) const{
	return allsuspects(
		get_xstart(bb.x), get_xend(bb.x),
		get_ystart(bb.y), get_yend(bb.y),
		get_zstart(bb.z), get_yend(bb.z));
}

vector<int> BoundingBox3DFinder::allsuspects(int x0, int x1, int y0, int y1, int z0, int z1) const{
	vector<int> ret;
	for (int ix=x0; ix<=x1; ++ix)
	for (int iy=y0; iy<=y1; ++iy)
	for (int iz=z0; iz<=z1; ++iz){
		addsuspects(ix, iy, iz, ret);
	}
	//no dublicates
	std::sort(ret.begin(), ret.end());
	auto iend = std::unique(ret.begin(), ret.end());
	ret.resize(iend - ret.begin());
	return ret;
}

void BoundingBox3DFinder::addsuspects(int ix, int iy, int iz, vector<int>& ret) const{
	int gi = iz*(mx+1)*(my+1)+iy*(mx+1)+ix;
	int sz = data[gi].size();
	ret.resize(ret.size()+sz);
	std::copy(data[gi].begin(), data[gi].end(), ret.end()-sz);
}
