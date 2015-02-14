#include "bgeom.h"
#include "addalgo.hpp"
#include "assert.h"

double Point::meas_section(const Point& p, const Point& L1, const Point& L2) noexcept{
	Vect a = L2-L1, b = p - L1;
	double ksi=vecDot(a,b)/vecDot(a,a);
	if (ksi>=1) return meas(p,L2);
	else if (ksi<=0) return meas(p,L1);
	else{
		double A[3] = { L1.y-L2.y, L2.x-L1.x, vecCrossZ(L1,L2) };
		double d0=A[0]*p.x+A[1]*p.y+A[2];
		d0*=d0; d0/=(A[0]*A[0]+A[1]*A[1]); //distance^2 to line 
		return d0;
	}
}

bool isOnSection(const Point& p, const Point& start, const Point& end, double& ksi, double eps){
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
		assert(!ISEQ(start.y, end.y));
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

vector<double> RefineSection(double a, double b, double Len, double Den){
	//std::cout<<"DUMMY REFINE SECTION"<<std::endl;
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





