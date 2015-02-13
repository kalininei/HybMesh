#include "contours.h"

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

void PContour::reverse_self(){
	std::reverse(pts.begin(), pts.end()); 
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

double PContour::area() const{
	double ret = 0;
	auto p1 = get_point(0);
	for (int i=1; i<n_points()-1; ++i){
		auto p2 = get_point(i);
		auto p3 = get_point(i+1);
		ret += triarea(*p1, *p2, *p3);
	}
	return ret;
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
	
int PContour::is_inside(const Point& p, const bool* inner_hint) const{
	//check whether contours is inner or outer
	bool is_inner = (inner_hint!=0) ? *inner_hint : (area()>0);
	//some random farpoint which lies outside any contour
	static Point farpoint = Point(1e3+1.23411+45.4*geps, 2.3e3+7.432+2.1*geps);
	//calculate number of crosses between [p, farpoint] and all contour segments
	int numcrosses = 0;
	double ksi[2];
	auto pc1=pts.back();
	for (auto pc2: pts){
		//if point lies on section stop loop
		if (isOnSection(p, *pc1, *pc2, ksi[0])) return 0;
		if (SectCross(p, farpoint, *pc1, *pc2, ksi)) ++numcrosses;
		pc1 = pc2;
	}
	//if even number of crosses -> inner point; else -> outer
	if ((is_inner && numcrosses % 2 == 1) || (!is_inner && numcrosses % 2 == 0)){
		return 1;
	} else {
		return -1;
	}
}

std::tuple<
	vector<Point*>,  //internal points
	vector<Point*>,  //points on contour
	vector<Point*>   //outer points
> PContour::filter_points(const vector<Point*>& points) const{
	//return value initialization
	std::tuple<vector<Point*>, vector<Point*>, vector<Point*>> ret;
	//check whether contours is inner or outer
	bool is_inner = (area()>0);
	//points loop
	for (auto p: points){
		int pos = is_inside(*p, &is_inner);
		if (pos==-1) std::get<2>(ret).push_back(p);
		if (pos==0) std::get<1>(ret).push_back(p);
		if (pos==1) std::get<0>(ret).push_back(p);
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

ContoursCollection::ContoursCollection(const vector<PContour>& cnts){
	for (auto c: cnts) add_contour(c);
}

void ContoursCollection::add_contour(const PContour& cnt){
	auto newc = aa::add_shared(contours, PContour(cnt));
	auto newe = aa::add_shared(entries, _entry(newc));
	auto p = cnt.get_point(0);
	auto tope = const_cast<ContoursCollection::_entry*>(efind(*p));
	if (tope!=NULL){
		//if new contour lie within existing one
		newe->upper = tope;
		tope->lower.push_back(newe);
	} else {
		//if new contour doesn't lie within any other
		//check if any top level contour lie within new one
		top_level.push_back(newe);
	}
	set_nesting();
}

const ContoursCollection::_entry* ContoursCollection::efind(const Point& p) const{
	for (auto e: top_level){
		auto r = e->find(p);
		if (r!=NULL) return r;
	}
	return NULL;
}

void ContoursCollection::set_nesting(){
	for (auto e: top_level) e->set_nesting(true);
}

int ContoursCollection::_entry::geom_inside(const Point& p) const{
	int r = data->is_inside(p, &is_inner);
	return (is_inner) ? r : -r;
}

void ContoursCollection::_entry::set_nesting(bool inner){
	if (is_inner != inner){
		is_inner = inner;
		data->reverse_self();
	}
	for (auto e: lower) e->set_nesting(!inner);
}

const ContoursCollection::_entry* ContoursCollection::_entry::find(const Point& p) const{
	if (geom_inside(p) == 1){
		for (auto e: lower){
			auto r1 = e->find(p);
			if (r1!=NULL) return r1;
		}
		return this;
	}
	return NULL;
}

vector<PContour> ContoursCollection::contours_list() const{
	vector<PContour> ret;
	for (auto sc: contours){
		ret.push_back(PContour(*sc));
	}
	return ret;
}

int ContoursCollection::is_inside(const Point& p) const{
	//point is within if it lies within any inner contour
	//and without all its child outer contours
	for (auto inner_e: entries) if (inner_e->is_inner){
		int r = inner_e->geom_inside(p);
		if (r==0) return 0;
		if (r==1){
			bool flag = true;
			for (auto outer_e: inner_e->lower){
				int r2=outer_e->geom_inside(p);
				if (r2==0) return 0;
				else if (r2==1) { flag = false; break; }
			}
			if (flag) return 1;
		}
	}
	return -1;
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
