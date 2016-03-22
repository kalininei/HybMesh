#include "contours.h"
#include <tuple>
#include "fileproc.h"
#include "assert.h"

PContour PContour::reverse() const {
	PContour ret(*this);
	std::reverse(ret.pts.begin(), ret.pts.end());
	return ret;
}

void PContour::reverse_self(){
	std::reverse(pts.begin(), pts.end());
}

PContour PContour::simplify() const {
	auto ret = PContour(*this);
	ret.simplify_self();
	return ret;
}

//simplify itself
void PContour::simplify_self(){
	std::set<int> simp_pnt;
	for (int i=0; i<n_points(); ++i) if (!is_corner_point(i)) simp_pnt.insert(i);
	aa::remove_entries(pts, simp_pnt);
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
	vector<const Point*> ret;
	bool is_inner = (area()>0);
	for (auto p: pts){
		int pos = is_inside(*p, &is_inner);
		if (pos!=-1) ret.push_back(p);
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

//find internal point: cross algorithm
Point PContour::inside_point_ca() const{
	CGBoundingBox bb(*this, 0.5);
	Point farp1(bb.xmin-134.11*geps, bb.ymin), farp2(bb.xmax+243.43*geps, bb.ymax);
	//get all crosses
	vector<double> crosses;
	double ksieta[] = {-1, -1};
	for (int i=0; i<n_points(); ++i){
		SectCross(farp1, farp2, *get_point(i), *get_point(i+1), ksieta);
		if (ksieta[1] > geps && ksieta[1] < 1+geps){
			crosses.push_back(ksieta[0]);
		}
	}
	if (crosses.size() < 2) throw std::runtime_error("failed to find inside point");
	//sort crosses
	sort(crosses.begin(), crosses.end());
	//find longest cross inside segment
	int ti = 0; double td = crosses[1] - crosses[0];
	for (int i=2; i<crosses.size()-1; i+=2){
		if (crosses[i+1]-crosses[i] > td){
			ti = i;
			td = crosses[i+1]-crosses[i];
		}
	}
	//return a point
	double w = (crosses[ti] + crosses[ti+1])/2.0;
	Point ret=Point::Weigh(farp1, farp2, w);
#ifndef NDEBUG
	bool hint = true;
	assert(this->is_inside(ret, &hint) == INSIDE);
#endif
	return ret;
}

vector<double> PContour::section_lenghts() const{
	vector<double> ret;
	for (int i=0; i<n_points(); ++i){
		ret.push_back(Point::dist(*get_point(i), *get_point(i+1)));
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

vector<double> PContour::max_chdist() const{
	vector<double> ret;
	vector<double> lens = section_lenghts();
	double dprev = lens.back();
	for (auto d: lens){
		ret.push_back(std::max(d,dprev));
		dprev = d;
	}
	return ret;
}

int PContour::is_inside(const Point& p, const bool* inner_hint) const{
	//check whether contours is inner or outer
	bool is_inner = (inner_hint!=0) ? *inner_hint : (area()>0);
	//simple bounding box check
	if (CGBoundingBox(*this).whereis(p) == OUTSIDE) return is_inner ? OUTSIDE : INSIDE;
	//some random farpoint which lies outside any contour
	static Point farpoint = Point(1e3+1.23411+45.4*geps, 2.3e3+7.432+2.1*geps);
	//calculate number of crosses between [p, farpoint] and all contour segments
	int numcrosses = 0;
	double ksi[2];
	auto pc1=pts.back();
	for (auto pc2: pts){
		//if point lies on section stop loop
		if (isOnSection(p, *pc1, *pc2, ksi[0])) return BOUND;
		if (SectCross(p, farpoint, *pc1, *pc2, ksi)) ++numcrosses;
		pc1 = pc2;
	}
	//if even number of crosses -> inner point; else -> outer
	if ((is_inner && numcrosses % 2 == 1) || (!is_inner && numcrosses % 2 == 0)){
		return INSIDE;
	} else {
		return OUTSIDE;
	}
}

bool PContour::is_corner_point(int i) const{
	auto p = get_point(i);
	auto pnext = get_point(i+1);
	auto pprev = get_point(i-1);
	double area3 = triarea(*pprev, *p, *pnext);
	return fabs(area3)>geps2;
}

bool PContour::is_corner_point(int i, double dev) const{
	auto p = get_point(i);
	auto pnext = get_point(i+1);
	auto pprev = get_point(i-1);
	auto a = Angle(*pprev, *p, *pnext)/M_PI*180;
	return fabs(a-180) > dev;
}

void PContour::delete_by_index(const std::set<int>& badind){
	aa::remove_entries(pts, badind);
}

std::tuple<
	vector<const Point*>,  //internal points
	vector<const Point*>,  //points on contour
	vector<const Point*>   //outer points
> PContour::filter_points(const vector<const Point*>& points) const{
	//return value initialization
	std::tuple<vector<const Point*>, vector<const Point*>, vector<const Point*>> ret;
	//check whether contours is inner or outer
	bool is_inner = (area()>0);
	//points loop
	for (auto p: points){
		int pos = is_inside(*p, &is_inner);
		if (pos==OUTSIDE) std::get<2>(ret).push_back(p);
		if (pos==BOUND) std::get<1>(ret).push_back(p);
		if (pos==INSIDE) std::get<0>(ret).push_back(p);
	}
	return ret;
}

vector<std::pair<double, Point>>
PContour::intersections(const PContour& p) const {
	vector<std::pair<double, Point>> ret;
	std::set<Point> pntset;
	CGBoundingBox bbox1(*this), bbox2(p);
	if (!bbox1.has_common_points(bbox2)) return ret;
	double ksieta[] = {-1,-1};
	for (int i=0; i<n_points(); ++i){
		auto a1=get_point(i), a2=get_point(i+1);
		if (!bbox2.contains(*a1, *a2)) continue;
		for (int j=0; j<p.n_points(); ++j){
			auto b1=p.get_point(j);
			auto b2=p.get_point(j+1);
			SectCross(*a1, *a2, *b1, *b2, ksieta);
			if (ksieta[0]>-geps && ksieta[0]<1+geps
					&& ksieta[1]>-geps && ksieta[1]<1+geps){
				auto newp =  Point::Weigh(*a1, *a2, ksieta[0]);
				if (pntset.insert(newp).second){
					ret.push_back(std::make_pair(i+ksieta[0],newp));
				}
			}

		}
	}
	return ret;
}

vector<std::pair<double, Point>>
PContour::intersections(const vector<PContour>& p) const{
	vector<std::pair<double, Point>> ret;
	for (int i=0; i<p.size(); ++i){
		auto ri = intersections(p[i]);
		for (auto r: ri) ret.push_back(r);
	}
	return ret;
}

ContoursCollection::ContoursCollection(const vector<PContour>& cnts){
	for (auto c: cnts) add_contour(c);
}

void ContoursCollection::add_contour(const PContour& cnt){
	auto newc = aa::add_shared(contours, PContour(cnt));
	auto newe = aa::add_shared(entries, _entry(newc));
	auto p = cnt.get_point(0);
	//find highest contour which contain new one
	auto tope = const_cast<ContoursCollection::_entry*>(efind(*p));
	if (tope!=NULL){
		//check if any tope contour lie within new one
		for (auto c: tope->lower){
			const Point* p = c->data->get_point(0);
			if (newe->geom_inside(*p)==INSIDE) newe->lower.push_back(c);
		}
		//delete newe internal contours from tope list
		for (auto c: newe->lower) tope->lower.remove(c);
		//add newe to tope as its lower contour
		tope->lower.push_back(newe);
	} else {
		//if new contour doesn't lie within any other
		//check if any top level contour lie within new one
		for (auto c: top_level){
			const Point* p = c->data->get_point(0);
			if (newe->geom_inside(*p)==INSIDE) newe->lower.push_back(c);
		}
		//delete newe internal contours from top level list
		for (auto c: newe->lower) top_level.remove(c);
		//add newe to top level
		top_level.push_back(newe);
	}

	//set upper for all entries
	for (auto e: entries){
		for (auto le: e->lower) le->upper = e.get();
	}
	for (auto e: top_level) e->upper = 0;
	//check directions
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
	if (geom_inside(p) == INSIDE){
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

int ContoursCollection::get_parent_index(int i) const{
	auto parent_entry = entries[i]->upper;
	if (parent_entry == 0) return -1;
	for (int i=0; i<entries.size(); ++i){
		if (parent_entry == entries[i].get()) return i;
	}
	throw std::runtime_error("Failed to find contours collection parent entry");
}

int ContoursCollection::is_inside(const Point& p) const{
	//point is within if it lies within any inner contour
	//and without all its child outer contours
	for (auto inner_e: entries) if (inner_e->is_inner){
		int r = inner_e->geom_inside(p);
		if (r==BOUND) return BOUND;
		if (r==INSIDE){
			bool flag = true;
			for (auto outer_e: inner_e->lower){
				int r2=outer_e->geom_inside(p);
				if (r2==0) return BOUND;
				else if (r2==INSIDE) { flag = false; break; }
			}
			if (flag) return INSIDE;
		}
	}
	return OUTSIDE;
}
std::tuple<
	vector<int>,  //internal points
	vector<int>,  //points on contour
	vector<int>   //outer points
> ContoursCollection::filter_points_i(const vector<Point>& points) const{
	vector<const Point*> p; p.reserve(points.size());
	for (auto& it: points) p.push_back(&it);
	return filter_points_i(p);
}

std::tuple<
	vector<int>,  //internal points
	vector<int>,  //points on contour
	vector<int>   //outer points
> ContoursCollection::filter_points_i(const vector<const Point*>& points) const{
	//return value
	std::tuple<vector<int>, vector<int>, vector<int>> ret;
	auto& ipnt = std::get<0>(ret);
	auto& cpnt = std::get<1>(ret);
	auto& opnt = std::get<2>(ret);
	//put simplifies contours into all entries
	simplify_entries();
	//loop over all points
	for (size_t i=0; i<points.size(); ++i){
		int r = is_inside(*points[i]);
		switch (r){
			case OUTSIDE: opnt.push_back(i); break;
			case INSIDE: ipnt.push_back(i); break;
			default: cpnt.push_back(i);
		}
	}
	//return original contours
	unsimplify_entries();

	return ret;
}

void ContoursCollection::simplify_entries() const{
	_origcont.clear();
	_simpcont.clear();
	for(auto& e: entries){
		//backup
		_origcont.push_back(e->data);
		//build new
		e->data = aa::add_shared(_simpcont, e->data->simplify());
		//check for adequate contour
		if (e->data->n_points() < 3){
			throw std::runtime_error("Invalid simplified contour");
		}
	}
}

void ContoursCollection::unsimplify_entries() const{
	auto it = _origcont.begin();
	for (auto& e: entries) e->data = *it++;
}

double ContoursCollection::area() const{
	auto clist = contours_list();
	double ret = 0;
	for (auto c: clist){
		ret += c.area();
	}
	return ret;
}

ContoursCollection ContoursCollection::level_01(int i) const{
	ContoursCollection res;
	res.add_contour(contour(i));
	for (auto cc: get_childs(i)) res.add_contour(*cc);
	return res;
}

ContoursCollection ContoursCollection::cut_by_level(int ist, int iend) const{
	ContoursCollection res;
	for (int i=0; i<n_cont(); ++i){
		int lv = get_level(i);
		if (lv>=ist && lv<=iend) res.add_contour(contour(i));
	}
	return res;
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

void PointsContoursCollection::do_scale(const ScaleBase& sc){
	sc.p_scale(pdata.begin(), pdata.end());
}
void PointsContoursCollection::undo_scale(const ScaleBase& sc){
	sc.p_unscale(pdata.begin(), pdata.end());
}

PointsContoursCollection::PointsContoursCollection(const vector<double>& pts, const vector<int>& eds){
	vector<Point> p;
	for (int i=0; i<pts.size(); i+=2){
		p.push_back(Point(pts[i], pts[i+1]));
	}
	build(p, eds);
}
PointsContoursCollection::PointsContoursCollection(const vector<Point>& pts, const vector<int>& eds){
	build(pts, eds);
}

void PointsContoursCollection::build(const vector<Point>& pts, const vector<int>& eds){
	for (auto p: pts) aa::add_shared(pdata, p);
	//filling edges
	for (int i=0; i<eds.size(); i+=2){
		edges.push_back(Edge(this, eds[i], eds[i+1], i/2));
	}
	//assemble edges into a closed contours
	std::list<Edge*> edset;
	for (auto& e: edges) edset.push_back(&e);
	vector<vector<Edge*>> clist;
	while (edset.size() > 0){
		clist.push_back(vector<Edge*>());
		auto& v = clist.back();
		auto pos0 = edset.begin();
		int cur = (**pos0).i1;
		int last = (**pos0).i0;
		v.push_back(*pos0);
		edset.pop_front();
		while (cur!=last){
			auto it = edset.begin();
			while (it != edset.end()){
				int i0 = (**it).i0;
				int i1 = (**it).i1;
				if (i0 == cur || i1 == cur){
					cur = (i0 == cur) ? i1 : i0;
					break;
				} else ++it;
			}
			if (it != edset.end()){
				v.push_back(*it);
				edset.erase(it);
			} else {
				throw std::runtime_error("failed to find closed contour");
			}
		}
		//edges renumbering
		for (int i=1; i<v.size(); ++i){
			if ((*v[i]).i0 != (*v[i-1]).i1){
				std::swap((*v[i]).i0, (*v[i]).i1);
			}
		}
	}
	//build PContours
	for (auto& c: clist){
		auto cont = PContour();
		for (auto& r: c){
			cont.add_point(pdata[(*r).i0].get());
		}
		ContoursCollection::add_contour(cont);
	}
}

PointsContoursCollection::PointsContoursCollection(const ContoursCollection& col){
	vector<Point> pts; vector<int> eds;
	for (int i=0; i<col.n_cont(); ++i){
		auto c = col.get_contour(i);
		int sz0 = pts.size();
		for (int j=0; j<c->n_points(); ++j){
			pts.push_back(*c->get_point(j));
			eds.push_back(pts.size()-1);
			eds.push_back(pts.size());
		}
		eds.back() = sz0;
	}
	build(pts, eds);
}

ContoursCollection PointsContoursCollection::shallow_copy(){
	ContoursCollection res;
	for (int i=0; i<n_cont(); ++i) res.add_contour(contour(i));
	return res;
}

vector<int> PointsContoursCollection::edge_correlation(
		const PointsContoursCollection& src,
		const PointsContoursCollection& tar){

	//group edges by their angles
	aa::DoubleMap<vector<const Edge*>> src_ang(geps);
	for (int i=0; i<src.n_edges(); ++i){
		double an = src.edges[i].get_angle();
		auto fnd = src_ang.emplace(an, vector<const Edge*>());
		fnd.first->second.push_back(&src.edges[i]);
	}

	//for central points of each target edges find an edge
	//in src among those with the same angle
	double ksi;
	vector<int> ret(tar.n_edges(), -1);
	for (int i=0; i<tar.n_edges(); ++i){
		auto& e = tar.edges[i];
		auto fnd = src_ang.find(e.get_angle());
		if (fnd!=src_ang.end()){
			Point tcp = e.center();
			for (auto& se: fnd->second){
				if (isOnSection(tcp, *se->pnt0(), *se->pnt1(), ksi)){
					ret[i] = se->index;
					break;
				}
			}
		}
	}
	return ret;
}

double PointsContoursCollection::Edge::get_angle() const{
	if (_angle < -0.5){
		Vect v = *pnt1() - *pnt0();
		_angle = atan2(v.y, v.x);
		if (_angle<0) _angle+=M_PI;
		if (ISEQ(_angle, M_PI)) _angle=0;
	}
	return _angle;
}

CGBoundingBox::CGBoundingBox(const vector<CGBoundingBox>& bb, double e):BoundingBox(){
	init();
	for (auto& b: bb){
		if (xmin > b.xmin) xmin = b.xmin;
		if (ymin > b.ymin) ymin = b.ymin;
		if (xmax < b.xmax) xmax = b.xmax;
		if (ymax < b.ymax) ymax = b.ymax;
	}
	widen(e);
}

CGBoundingBox::CGBoundingBox(const PContour& cont, double e):BoundingBox(){
	init();
	for (int i=0; i<cont.n_points(); ++i) WidenWithPoint(*cont.get_point(i));
	widen(e);
}

CGBoundingBox::CGBoundingBox(const ContoursCollection& col, double e):BoundingBox(){
	init();
	for (int i=0; i<col.n_cont(); ++i){
		auto c = col.get_contour(i);
		for (int j=0; j<c->n_points(); ++j) WidenWithPoint(*c->get_point(j));
	}
	widen(e);
}

Contour CGBoundingBox::get_contour() const {
	vector<Point> p {Point(xmin, ymin),
		Point(xmax, ymin), Point(xmax, ymax), Point(xmin, ymax)};
	return Contour(p);
}


CGBoundingBox::CGBoundingBox(const Point& p1, const Point& p2, double e):BoundingBox(){
	xmin = std::min(p1.x, p2.x);
	ymin = std::min(p1.y, p2.y);
	xmax = std::max(p1.x, p2.x);
	ymax = std::max(p1.y, p2.y);
	widen(e);
}
