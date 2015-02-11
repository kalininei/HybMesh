#include "bgeom.h"
#include "addalgo.hpp"

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

void PContour::select_points(const std::vector<Point*>& pts, 
		std::vector<Point*>& inner, std::vector<Point*>& outer) const{
	auto Rect = rectangle_bnd();
	for (auto p: pts){
		//check if point is within rectangle
		//if (p->x < Rect.first.x ||
				//p->x > Rect.second.x ||
				//p->y < Rect.first.y ||
				//p->y > Rect.second.y){
		if (ISLOWER(p->x, Rect.first.x) ||
				ISGREATER(p->x, Rect.second.x) ||
				ISLOWER(p->y, Rect.first.y) ||
				ISGREATER(p->y, Rect.second.y)){
			outer.push_back(p);
		}else{
			//TODO
			inner.push_back(p);
		}
	}
}
PContour PContour::reverse() const { 
	PContour ret(*this);
	std::reverse(ret.pts.begin(), ret.pts.end()); 
	return ret;
}

vector<double> PContour::meas_points(const vector<const Point*>& pts) const {
	auto inner_pts_vec = find_inner(pts);
	auto inner_pts = std::set<const Point*>(inner_pts_vec.begin(), inner_pts_vec.end());
	vector<double> ret;
	for(auto p: pts){
		double v = meas_to_point(*p);
		if (inner_pts.find(p)==inner_pts.end()) v*=-1;
		ret.push_back(v);
	}
	return ret;
}

//double PContour::extent() const{
	//double ret = 0;
	//for (int i=0; i<n_points(); ++i){
		//for (int j=i+1; j<n_points(); ++j){
			//double d = Point::meas(*pts[i], *pts[j]);
			//if (d>ret) ret = d;
		//}
	//}
	//return sqrt(ret);
//}

double PContour::meas_to_point(const Point& p) const{
	auto pprev = pts.back();
	double ret = 1e100;
	for (auto pc: pts){
		double v = Point::meas_section(p, *pprev, *pc);
		if (v<ret) ret = v;
		pprev = pc;
	}
	return ret;
}

vector<const Point*> PContour::find_inner(const vector<const Point*>& pts) const{
	std::cout<<"Dummy find_inner"<<std::endl;
	vector<const Point*> ret;
	for (auto p: pts){
		if (ISEQGREATER(p->x, 0.6) || ISEQGREATER(p->y, 0.6) ||
				ISEQLOWER(p->x, 0.3) || ISEQLOWER(p->y, 0.3)){}
		else{
			ret.push_back(p);
		}
	}
	return ret;
}
//vector<const Point*> PContour::filter_points(const vector<const Point*>& pts, double dist1, double dist2) const{
	//vector<const Point*> inner_pts=find_inner(pts);
	//if (dist1<=0 && dist2 > extent()) return inner_pts;
	//vector<const Point*> ret;
	//dist1*=dist1; dist2*=dist2;
	//for (auto p: inner_pts){
		//double m = meas_to_point(*p);
		//if (m>=dist1 && m<=dist2) ret.push_back(p);
	//}
	//return ret;
//}

Contour PContour::widen_contour(double buffer_size) const{
	//TODO
	return Contour(*this);
}

vector<double> PContour::section_lenghts() const{
	vector<double> ret;
	auto pprev = pts.back();
	for (auto p: pts){
		ret.push_back(Point::dist(*pprev,*p));
		pprev = p;
	}
	return ret;
}
vector<double> PContour::chdist() const{
	vector<double> ret;
	vector<double> lens = section_lenghts();
	double dprev = lens.back();
	for (auto d: lens){
		ret.push_back((d+dprev)/2.0);
		dprev = d;
	}
	return ret;
}

std::pair<Point, Point> PContour::rectangle_bnd() const{
	double minx = pts[0]->x, maxx = pts[0]->x;
	double miny = pts[0]->y, maxy = pts[0]->y;
	for (auto p: pts){
		if (p->x<minx) minx = p->x;
		if (p->y<miny) miny = p->y;
		if (p->x>maxx) maxx = p->x;
		if (p->y>maxy) maxy = p->y;
	}
	return std::make_pair(Point(minx, miny), Point(maxx, maxy));
	
}

Contour::Contour(const PContour& c): PContour(){
	for (int i=0; i<c.n_points(); ++i){
		add_point(*c.get_point(i));
	}
}

Contour::Contour(const std::vector<Point>& _pts): PContour(){
	for (auto& p: _pts) add_point(p);
}

void Contour::add_point(Point* p){
	Point pcopy(*p);
	add_point(pcopy);
}

void Contour::add_point(const Point& p){
	auto pp = aa::add_shared(pdata, p);
	PContour::add_point(pp);
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





