#include "hmgrid3dgeom.hpp"
using namespace HMGrid3D;

BoundingBox3D::BoundingBox3D(const Cell& c){
	auto v = c.allvertices();
	fill(v);
}

BoundingBox3D::BoundingBox3D(const Surface& c){
	auto v = c.allvertices();
	fill(v);
}
BoundingBox3D::BoundingBox3D(const VertexData& v){
	fill(v);
}

void BoundingBox3D::fill(const ShpVector<Vertex>& vv){
	if (vv.size() == 0){
		xmin = xmax = ymin = ymax = zmin = zmax = 0;
	} else {
		xmin = xmax = vv[0]->x;
		ymin = ymax = vv[0]->y;
		zmin = zmax = vv[0]->z;
		for (int i=1; i<vv.size(); ++i){
			if (vv[i]->x > xmax) xmax = vv[i]->x;
			if (vv[i]->x < xmin) xmin = vv[i]->x;
			if (vv[i]->y > ymax) ymax = vv[i]->y;
			if (vv[i]->y < ymin) ymin = vv[i]->y;
			if (vv[i]->z > zmax) zmax = vv[i]->z;
			if (vv[i]->z < zmin) zmin = vv[i]->z;
		}
	}
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
	int y1 = pos1D(bb.ymin, ymin, ymax);
	int y2 = pos1D(bb.ymax, ymin, ymax);
	int z1 = pos1D(bb.zmin, zmin, zmax);
	int z2 = pos1D(bb.zmax, zmin, zmax);

	if (!x1 && !x2 && !y1 && !y2 && !z1 && !z2) return INSIDE;
	else if ((x1 == x2 && x1 != 0) || (y1 == y2 && y1 != 0) || (z1 == z2 && z1 != 0)) return OUTSIDE;
	else return BOUND;
}
int BoundingBox3D::position(const Vertex& v) const{
	if (ISLOWER(v.x, xmin) || ISGREATER(v.x, xmax) ||
		ISLOWER(v.y, ymin) || ISGREATER(v.y, ymax) ||
		ISLOWER(v.z, zmin) || ISGREATER(v.z, zmax)) return OUTSIDE;
	else if (ISGREATER(v.x, xmin) && ISLOWER(v.x, xmax) &&
		ISGREATER(v.y, ymin) && ISLOWER(v.y, ymax) &&
		ISGREATER(v.z, zmin) && ISLOWER(v.z, zmax)) return INSIDE;
	return BOUND;
}

//calculate measures
double BoundingBox3D::meas(const BoundingBox3D& bb) const{
	auto di = [](double x, double y, double a, double b)->double{
		int p1 = pos1D(x, a, b);
		int p2 = pos1D(y, a, b);
		if (p1 < 0 && p2 < 0) return a-y;
		else if (p1 > 0 && p2 > 0) return x-b;
		else return 0.0;
	};

	double dx = di(xmin, xmax, bb.xmin, bb.xmax);
	double dy = di(ymin, ymax, bb.ymin, bb.ymax);
	double dz = di(zmin, zmax, bb.zmin, bb.zmax);

	return dx*dx + dy*dy + dz*dz;
}

//edit
BoundingBox3D BoundingBox3D::widen(double eps) const{
	BoundingBox3D ret(*this);
	ret.xmin-=eps; ret.xmax+=eps;
	ret.ymin-=eps; ret.ymax+=eps;
	ret.zmin-=eps; ret.zmax+=eps;
	return ret;
}

BoundingSphere3D  BoundingBox3D::inscribe_sphere() const{
	BoundingSphere3D ret;
	ret.x = (xmax + xmin)/2.;
	ret.y = (ymax + ymin)/2.;
	ret.z = (zmax + zmin)/2.;
	ret.rad = std::min(xmax-xmin, std::min(ymax-ymin, zmax-zmin));
	ret.rad/=2.0;
	return ret;
}


BoundingSphere3D::BoundingSphere3D(const Surface& c, BoundingBox3D* helper){
	auto v = c.allvertices();
	fill(v, helper);
}

int BoundingSphere3D::position(const BoundingBox3D& bb) const {
	auto di = [](double x, double a, double b)->double{
		int p1 = pos1D(x, a, b);
		if (p1 < 0) return a-x;
		else if (p1 > 0) return x-b;
		else return std::min(b-x, x-a);
	};
	int cpos = bb.position(Vertex(x, y, z));
	if (cpos == BOUND) return BOUND; 
	double dx = di(x, bb.xmin, bb.xmax);
	double dy = di(y, bb.ymin, bb.ymax);
	double dz = di(z, bb.zmin, bb.zmax);
	if (cpos == OUTSIDE){
		double dist2 = dx*dx + dy*dy + dz*dz;
		if (ISGREATER(dist2, rad*rad)) return OUTSIDE;
		else return BOUND;
	} else {
		double dist = std::min(dx, std::min(dy, dz));
		if (ISLOWER(rad, dist)) return INSIDE;
		else return BOUND;
	}
}

int BoundingSphere3D::position(const Vertex& v) const{
	double m = cmeas(v);
	if (ISGREATER(m, rad*rad)) return OUTSIDE;
	else if (ISLOWER(m, rad*rad)) return INSIDE;
	else return BOUND;
}

double BoundingSphere3D::meas(const BoundingBox3D& bb) const{
	auto di = [](double x, double a, double b)->double{
		int p1 = pos1D(x, a, b);
		if (p1 < 0) return a-x;
		else if (p1 > 0) return x-b;
		else return 0;
	};
	int cpos = bb.position(Vertex(x, y, z));
	if (cpos == OUTSIDE){
		double dx = di(x, bb.xmin, bb.xmax);
		double dy = di(y, bb.ymin, bb.ymax);
		double dz = di(z, bb.zmin, bb.zmax);
		double ret = sqrt(dx*dx + dy*dy + dz*dz) - rad;
		return ret*ret;
	} else return 0.0;
}

double BoundingSphere3D::cmeas(const Vertex& v) const{
	return ((v.x-x)*(v.x-x) + (v.y-y)*(v.y-y) + (v.z-z)*(v.z-z));
}

BoundingSphere3D BoundingSphere3D::widen(double eps) const{
	BoundingSphere3D ret(*this);
	ret.rad += eps;
	return ret;
}
void BoundingSphere3D::fill(const ShpVector<Vertex>& v, BoundingBox3D* helper){
	//crude approximation
	BoundingBox3D* bb;
	if (helper == 0) bb = new BoundingBox3D(v);
	else bb = helper;
	auto bs = bb->inscribe_sphere();
	for (auto& p: v){
		if (bs.position(*p) == OUTSIDE){
			double rnew = (bs.x-p->x)*(bs.x-p->x);
			rnew += (bs.y-p->y)*(bs.y-p->y);
			rnew += (bs.z-p->z)*(bs.z-p->z);
			bs.rad = sqrt(rnew);
		}
	}
	rad = bs.rad;
	x = bs.x; y = bs.y; z = bs.z;
	if (helper == 0) delete bb;
}
