#include "bgeom2d.h"
#include "addalgo.hpp"
#include "assert.h"

double Point::meas_section(const Point& p, const Point& L1, const Point& L2, double& ksi) noexcept{
	Vect a = L2-L1, b = p - L1;
	if (L1==L2) return Point::meas(p, L1);
	ksi=vecDot(a,b)/vecDot(a,a);
	if (ksi>=1) { ksi=1.0; return meas(p,L2);}
	else if (ksi<=0) { ksi=0.0; return meas(p,L1);}
	else{
		double A[3] = { L1.y-L2.y, L2.x-L1.x, vecCrossZ(L1,L2) };
		double d0=A[0]*p.x+A[1]*p.y+A[2];
		d0*=d0; d0/=(A[0]*A[0]+A[1]*A[1]); //distance^2 to line 
		return d0;
	}
}

double Point::meas_section(const Point& p, const Point& L1, const Point& L2) noexcept{
	static double k;
	return meas_section(p, L1, L2, k);
}

bool isOnSection(const Point& p, const Point& start, const Point& end, double& ksi, double eps) noexcept{
	ksi = gbig;
	//check if p is ouside section square
	if ( (p.x+eps < start.x && p.x+eps < end.x) ||
	     (p.x-eps > start.x && p.x-eps > end.x) ||
	     (p.y+eps < start.y && p.y+eps < end.y) ||
	     (p.y-eps > start.y && p.y-eps > end.y) ) return false;
	//calculation
	Vect v1=(end-start);
	vecNormalize(v1);
	double Dor=v1.y*start.x-v1.x*start.y;
	double Dp =v1.y*p.x-v1.x*p.y;
	if (fabs(Dor-Dp)>eps) return false;
	Point pproj;

	//exact match
	if (fabs(Dor-Dp)<geps) pproj=p;
	else {
	//~eps match
		Vect v2(-v1.y, v1.x);
		v2*=(Dor-Dp); //since |v2|=|v1|=1.0
		pproj=p-v2;
	}

	if (ISEQ(start.x, end.x)){
		if(ISEQ(start.y, end.y)) {ksi=0.5; return p==start;}
		ksi=(pproj.y-start.y)/(end.y-start.y);
	}else{
		ksi=(pproj.x-start.x)/(end.x-start.x);
	}
	return (ksi>-geps && ksi<1+geps);
}

double determinant_2x2(double* A){
	return A[0]*A[3]-A[1]*A[2];
}

bool mat_inverse_2x2(double* A, double* B){
	double det=determinant_2x2(A);
	if (ISZERO(det)) return false;
	B[0]= A[3]/det;
	B[1]=-A[1]/det;
	B[2]=-A[2]/det;
	B[3]= A[0]/det;
	return true;
}

bool SectCross(const Point& p1S, const Point& p1E, const Point& p2S, const Point& p2E, double* ksieta) noexcept{
	double A[4]={p1E.x-p1S.x,  p2S.x-p2E.x,
		  p1E.y-p1S.y,  p2S.y-p2E.y};
	double B[4];
	if (!mat_inverse_2x2(A,B)){ 
		ksieta[0]=gbig; 
		ksieta[1]=gbig;
	} else {
		double b[2]={p2S.x-p1S.x,
			  p2S.y-p1S.y};
		ksieta[0]=B[0]*b[0]+B[1]*b[1];
		ksieta[1]=B[2]*b[0]+B[3]*b[1];
	}
	return (ksieta[0]>=0 && ksieta[0]<=1 && ksieta[1]>=0 && ksieta[1]<=1); 
}

bool SectCrossWRenorm(const Point& p1S, const Point& p1E, const Point& p2S, const Point& p2E, double* ksieta) noexcept{
	vector<Point> p {p1S, p1E, p2S, p2E};
	ScaleBase s = ScaleBase::doscale(p);
	bool ret = SectCross(p[0], p[1], p[2], p[3], ksieta);
	if (ksieta[0]<gbig && ksieta[1]<gbig){
		ksieta[0] *= s.L;
		ksieta[1] *= s.L;
	}
	return ret;
}

int LinePointWhereIs(const Point& p, const Point& L1, const Point& L2) noexcept{
	if (L1 == L2) return -1;
	double A[3] = { L1.y-L2.y, L2.x-L1.x, vecCrossZ(L1,L2) };
	double v = A[0]*p.x + A[1]*p.y + A[2];
	if (ISZERO(v)) return 1;
	else return (v>0) ? 0 : 2;
}

vector<double> RefineSection(double a, double b, double Len, double Den){
	//get number of points
	double amean = (a+b)/2.0;
	double hmean = 2*a*b/(a+b);
	double gmean = sqrt(a*b);
	double h = (1-Den)*amean+Den*std::min(hmean, gmean);
	h = Len/std::ceil(Len/h);
	//refinement algo
	vector<double> ret;
	for (double s = h; s<Len-geps; s+=h) ret.push_back(s);
	return ret;
}



NodeFinder::NodeFinder(Point p0, double _lx, double _ly, int _nx, int _ny, double _eps):
		Nx(_nx), Ny(_ny), eps(1.1*_eps), 
		x0(p0.x-eps), y0(p0.y-eps), 
		x1(x0+_lx+2*eps), y1(y0+_ly+2*eps), 
		hx((x1-x0)/Nx), hy((y1-y0)/Ny), data(Nx*Ny){}
NodeFinder::NodeFinder(const std::pair<Point, Point>& rect, int _nx, int _ny, double _eps):
		Nx(_nx), Ny(_ny), eps(1.1*_eps), 
		x0(rect.first.x-eps), y0(rect.first.y-eps), 
		x1(rect.second.x+eps), y1(rect.second.y+2*eps), 
		hx((x1-x0)/Nx), hy((y1-y0)/Ny), data(Nx*Ny){}


bool NodeFinder::is_equal_point(const Point* p1, const Point* p2) const{
	return (fabs(p1->x-p2->x)<eps && (fabs(p1->y-p2->y)<eps));
}

const Point* NodeFinder::add(const Point* p){
	auto ind = get_index(p);
	if (ind.size()==0) return 0;
	//try to find point
	for (auto i: ind){
		for (auto cand: data[i]){
			if (is_equal_point(p, cand)) return cand;
		}
	}
	//add point
	for (auto i: ind) data[i].push_back(p);
	return p;
}

const Point* NodeFinder::find(const Point* p) const{
	auto ind = get_index(p);
	for (auto i: ind){
		for (auto cand: data[i]){
			if (is_equal_point(p, cand)) return cand;
		}
	}
	return 0;
}

vector<int> NodeFinder::get_index(const Point* p) const{
	if (p->x < x0 || p->y < y0 || p->x > x1 || p->y > y1) return vector<int>();
	double dx = p->x - x0, dy = p->y - y0;
	auto pure_find = [&](double x, double y)->int{
		int i = (int)std::floor(x/this->hx), j = (int)std::floor(y/this->hy);
		return (i>=0 && i<this->Nx && j>=0 && j<this->Ny) ? 
			this->to_glob(i, j) :
			-1;
	};
	int i1 = pure_find(dx-eps, dy-eps);
	int i2 = pure_find(dx-eps, dy+eps);
	int i3 = pure_find(dx+eps, dy+eps);
	int i4 = pure_find(dx+eps, dy-eps);
	vector<int> ret; ret.reserve(4);
	if (i1!=-1) ret.push_back(i1);
	if (i2!=-1 && i2!=i1) ret.push_back(i2);
	if (i3!=-1 && i3!=i2 && i3!=i1) ret.push_back(i3);
	if (i4!=-1 && i4!=i3 && i4!=i2 && i4!=i1) ret.push_back(i4);
	return ret;
}

void BoundingBox::init(){
	xmin =  gbig; xmax = -gbig;
	ymin =  gbig; ymax = -gbig;
}

void BoundingBox::widen(double e){
	xmin -= e; ymin -= e;
	xmax += e; ymax += e;
}

void BoundingBox::WidenWithPoint(const Point& p){
	if (p.x < xmin) xmin = p.x;
	if (p.x > xmax) xmax = p.x;
	if (p.y < ymin) ymin = p.y;
	if (p.y > ymax) ymax = p.y;
}

Point BoundingBox::Center() const{
	return Point((xmin+xmax)/2, (ymin+ymax)/2);
}

BoundingBox::BoundingBox(const vector<BoundingBox>& bb, double e){
	init();
	for (auto& b: bb){
		if (xmin > b.xmin) xmin = b.xmin;
		if (ymin > b.ymin) ymin = b.ymin;
		if (xmax < b.xmax) xmax = b.xmax;
		if (ymax < b.ymax) ymax = b.ymax;
	}
	widen(e);
}

BoundingBox::BoundingBox(const Point& p1, const Point& p2, double e){
	xmin = std::min(p1.x, p2.x);
	ymin = std::min(p1.y, p2.y);
	xmax = std::max(p1.x, p2.x);
	ymax = std::max(p1.y, p2.y);
	widen(e);
}


BoundingBox::BoundingBox(const Point& p1, double e){
	xmin = p1.x; xmax = p1.x;
	ymin = p1.y; ymax = p1.y;
	widen(e);
}


double BoundingBox::area() const{
	return (xmax-xmin)*(ymax-ymin);
}

int BoundingBox::whereis(const Point& p) const{
	if (ISLOWER(p.x, xmax) && ISLOWER(p.y, ymax) &&
			ISGREATER(p.x, xmin) && ISGREATER(p.y, ymin)) return INSIDE;
	if (ISGREATER(p.x, xmax) || ISGREATER(p.y, ymax) ||
			ISLOWER(p.x, xmin) || ISLOWER(p.y, ymin)) return OUTSIDE;
	return BOUND;
}

bool BoundingBox::has_common_points(const BoundingBox& bb) const{
	auto res=BoundingBox({*this, bb});
	auto sum=BoundingBox(0,0,lenx()+bb.lenx(), leny()+bb.leny());
	//if one lies outside the other
	if (ISGREATER(res.lenx(), sum.lenx()) ||
		ISGREATER(res.leny(), sum.leny())) return false;

	//if one lies inside the other
	//if (ISEQ(lenx(), res.lenx()) && ISEQ(leny(), res.lenx())){
	//        return whereis(Point(bb.xmin, bb.ymin)) == INSIDE &&
	//                whereis(Point(bb.xmax, bb.ymax)) == INSIDE;
	//}
	//if (ISEQ(bb.lenx(), res.lenx()) && ISEQ(bb.leny(), res.leny())){
	//        return bb.whereis(Point(xmin, ymin)) == INSIDE &&
	//                bb.whereis(Point(xmax, ymax)) == INSIDE;
	//}
	//has intersections
	return true;
}

bool BoundingBox::contains(const Point& p1, const Point& p2) const{
	if (whereis(p1) == INSIDE || whereis(p2) == INSIDE) return true;
	if (has_common_points(BoundingBox(p1,p2)) == false) return false;

	double A=p2.y-p1.y, B=p2.x-p1.x, C=p1.y*p2.x-p2.y*p1.x;
	auto func = [&](double x, double y)->int{ return SIGN(A*x+B*y+C); };

	int s = func(xmin, ymax);
	return !(s == func(xmax, ymin) && s == func(xmax, ymax) && s == func(xmin, ymax));

}

vector<Point> BoundingBox::FourPoints() const{
	return { Point(xmin, ymin), Point(xmax, ymin),
		Point(xmax, ymax), Point(xmin, ymax) };
}
