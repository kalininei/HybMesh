#include "cont_repart.hpp"
using namespace HMCont2D;

namespace{

struct data{
	vector<Point*> cont;
	vector<double> clen;
	std::list<double> div;
	double h;
	double hdiamond;
	double snap;

	//methods
	std::map<double, int> _finder;
	void build_finder(){
		if (_finder.size()==0){
			for (int i=0; i<clen.size(); ++i) _finder.emplace(clen[i], i);
		}
	}
	bool is_closed(){ return cont[0] == cont.back(); }
	std::list<double>::iterator insert_half(
			std::list<double>::iterator i0,
			std::list<double>::iterator i1){
		if (*i1 == 0) i1 = std::next(i0); //for closed
		return div.insert(i1, (*i0 + *i1)/2.0);
	}
	std::list<double>::iterator insert_half_wsnap(
			std::list<double>::iterator i0,
			std::list<double>::iterator i1){
		build_finder();
		if (*i1 == 0) i1 = std::next(i0); //for closed
		double val = (*i0 + *i1)/2.0;
		auto fnd = _finder.lower_bound(val);
		double vforw = fnd->first;
		double vback = (fnd != _finder.begin()) ? std::prev(fnd)->first : -1;
		if (vforw <= *i0 || vforw >= *i1) vforw = -1;
		else if (fabs(val - vforw) > snap*h) vforw = -1;
		if (vback >=0){
			if (vback <= *i0 || vback >= *i1) vback = -1;
			else if (fabs(val - vback) > snap*h) vback = -1;
		}
		if (vback<0 && vforw <0) return div.insert(i1, val);
		else if (vforw >=0 && vback <=0) return div.insert(i1, vforw);
		else if (vback >=0 && vforw <=0) return div.insert(i1, vback);
		else if (vforw < vback) return div.insert(i1, vforw);
		else return div.insert(i1, vback);
	}
	double dist(double w1, double w2){
		if (w2 == 0) w2 = clen.back();
		return w2 - w1;
	}
	double forw_length(std::list<double>::iterator it){
		if (++it != div.end()) return *it;
		else{
			assert(is_closed());
			return *div.begin();
		}
	}
	double backw_length(std::list<double>::iterator it){
		if (it != div.begin()) return *(--it);
		else{
			assert(is_closed());
			return *(----div.end());
		}
	}

	Point get_point(double w1){
		build_finder();
		auto f1 = _finder.lower_bound(w1);
		assert(f1 != _finder.end());
		if (ISEQ(clen[f1->second], w1)){
			return *cont[f1->second];
		} else {
			int fstart = f1->second;
			double l1 = clen[fstart-1];
			double l2 = clen[fstart];
			double t = (w1-l1)/(l2-l1);
			return Point::Weigh(*cont[fstart-1], *cont[fstart], t);
		}
	}

	vector<Point> get_segment(double w1, double w2){
		if (w2 == 0) w2 = clen.back();
		assert(w2>w1);
		build_finder();
		vector<Point> ret;
		//first point
		auto f1 = _finder.lower_bound(w1);
		assert(f1 != _finder.end());
		int fstart;
		if (ISEQ(clen[f1->second], w1)){
			fstart = f1->second + 1;
			ret.push_back(*cont[fstart-1]);
		} else {
			fstart = f1->second;
			double l1 = clen[fstart-1];
			double l2 = clen[fstart];
			double t = (w1-l1)/(l2-l1);
			ret.push_back(Point::Weigh(*cont[fstart-1], *cont[fstart], t));
		}
		//last point
		auto f2 = _finder.lower_bound(w2);
		assert(f2 != _finder.end());
		int fend;
		Point rback;
		if (ISEQ(clen[f2->second], w2)){
			fend = f2->second-1;
			rback = *cont[fend+1];
		} else {
			fend = f2->second-1;
			double l1 = clen[fend];
			double l2 = clen[fend+1];
			double t = (w2-l1)/(l2-l1);
			rback = Point::Weigh(*cont[fend], *cont[fend+1], t);
		}
		for (int i=fstart; i<=fend; ++i){
			ret.push_back(*cont[i]);
		}
		ret.push_back(rback);
		return ret;
	}
};

void add_angle_points(data& d, double zeroangle){
	for (int i=1; i<d.clen.size()-1; ++i){
		int iprev = i-1, inext = i+1;
		inext = i+1;
		double a = Angle(*d.cont[i-1], *d.cont[i], *d.cont[i+1]);
		if (ISEQLOWER(a, M_PI - zeroangle) || ISEQGREATER(a, M_PI + zeroangle)){
			d.div.push_back(d.clen[i]);
		}
	}
	d.div.sort();
	d.div.unique();
}

void uniform_division(data& d){
	auto it0 = d.div.begin();
	auto it1 = std::next(it0);
	while(it1 != d.div.end()){
		double w = *it1 - *it0;
		int n = std::max(1.0, round(w/d.h));
		double hn = w/n;
		for (int i=1; i<n; ++i) d.div.insert(it1, *it0 + i*hn);
		it0 = it1++;
	}
}

void sortout_points_within_triangle(const Point& p1, const Point& p2, const Point& p3,
		std::set<Point*>& points){
	double j11 = p2.x - p1.x, j21 = p3.x - p1.x;
	double j12 = p2.y - p1.y, j22 = p3.y - p1.y;
	double modj = (j22*j11 - j21*j12);

	if (fabs(modj) < geps*geps){
		double ksi;
		auto it = points.begin();
		while (it!=points.end()){
			if (isOnSection(**it, p1, p2, ksi) ||
			    isOnSection(**it, p2, p3, ksi) ||
			    isOnSection(**it, p3, p1, ksi)){
				it = points.erase(it);
			} else ++it;
		}
	} else {
		double ksi, eta;
		auto it = points.begin();
		while (it!=points.end()){
			auto p = *it;
			ksi = ( j22*(p->x - p1.x) - j21*(p->y - p1.y))/modj;
			eta = (-j12*(p->x - p1.x) + j11*(p->y - p1.y))/modj;
			if (ksi<-geps || ksi>1+geps || eta<-geps || eta>1-ksi+geps) ++it;
			else it = points.erase(it);
		}
	}
}

bool check_diamond(data& d, double w1, double w2){
	vector<Point> seg = d.get_segment(w1, w2);
	if (seg.size() == 2) return true;
	std::set<Point*> pset;
	for (int i=1; i<seg.size()-1; ++i) pset.insert(&seg[i]);
	Vect pv = seg.back() - seg[0];
	pv.set(-pv.y*d.hdiamond, pv.x*d.hdiamond);
	Point p3 = (seg[0] + seg.back())/2.0 + pv;
	sortout_points_within_triangle(seg[0], seg.back(), p3, pset);
	if (pset.size() == 0) return true;
	p3 = (seg[0] + seg.back())/2.0 - pv;
	sortout_points_within_triangle(seg[0], seg.back(), p3, pset);
	return (pset.size()==0);
}

void diamond(data& d, std::list<double>::iterator it1, std::list<double>::iterator it2){
	if (!check_diamond(d, *it1, *it2)){
		auto it = d.insert_half_wsnap(it1, it2);
		diamond(d, it1, it);
		diamond(d, it, it2);
	}
}

void diamond_rule(data& d){
	auto it0 = d.div.begin();
	auto it1 = std::next(it0);
	while (it1 != d.div.end()){
		diamond(d, it0, it1);
		it0 = it1++;
	}
}

void smooth_length(data& d,
		std::list<double>::iterator it0,
		std::list<double>::iterator it1,
		std::list<double>::iterator it2){
	double l1 = d.dist(*it0, *it1);
	double l2 = d.dist(*it1, *it2);
	if (l2 > 3*l1){
		auto it = d.insert_half_wsnap(it1, it2);
		smooth_length(d, it0, it1, it);
	} else if (l1 > 3*l2){
		auto it = d.insert_half_wsnap(it0, it1);
		smooth_length(d, it, it1, it2);
	}
}

void one_third_rule(data& d){
	if (d.div.size() <= 2) return;
	auto it0 = (d.is_closed()) ? ----d.div.end() : d.div.begin();
	auto it1 = (d.is_closed()) ? d.div.begin() : std::next(it0);
	auto it2 = std::next(it1);
	while (it2 != d.div.end()){
		smooth_length(d, it0, it1, it2);
		it1++;
		it2 = std::next(it1);
		it0 = std::prev(it1);
	}
}

vector<vector<int>> build_distrib(data& d){
	d.build_finder();
	vector<vector<int>> ret;
	auto it0 = d.div.begin();
	auto it1 = std::next(it0);
	while(it1 != d.div.end()){
		ret.push_back(vector<int>());
		auto& r = ret.back();
		auto fnd1 = d._finder.lower_bound(*it0);
		int fstart;
		if (ISEQ(fnd1->first, *it0)){
			r.push_back(fnd1->second);
			fstart = fnd1->second + 1;
		} else {
			r.push_back(-1);
			fstart = fnd1->second;
		}
		auto fnd2 = d._finder.lower_bound(*it1);
		int fend;
		double rback;
		if (ISEQ(fnd2->first, *it1)){
			rback = fnd2->second;
			fend = fnd2->second - 1;
		} else{
			rback = -1;
			fend = fnd2->second - 1;
		}
		for (int i=fstart; i<=fend; ++i) r.push_back(i);
		r.push_back(rback);
		it0 = it1++;
	}
	return ret;
}

bool try_forward_move(data& d, std::list<double>::iterator it,
		vector<int>& b, vector<int>& f){
	assert(f[0] == -1 && b.back() == -1 && f.size()>2);
	double w1 = d.backw_length(it);
	double w2 = d.clen[f[1]];
	if (check_diamond(d, w1, w2)){
		b.back() = f[1];
		*it = d.clen[f[1]];
		f.erase(f.begin());
		return true;
	} else {
		return false;
	}
}
bool try_backward_move(data& d, std::list<double>::iterator it,
		vector<int>& b, vector<int>& f){
	assert(f[0] == -1 && b.back() == -1 && b.size()>2);
	double w1 = d.clen[b[b.size()-2]];
	double w2 = d.forw_length(it);
	if (check_diamond(d, w1, w2)){
		b.resize(b.size()-1);
		*it = d.clen[b.back()];
		f[0] = b.back();
		return true;
	} else {
		return false;
	}
}

vector<vector<int>> snapping(data& d){
	vector<vector<int>> ret = build_distrib(d);
	auto it1 = d.div.begin();
	auto retf = ret.begin();
	for (int i=0; i<ret.size(); ++i, ++it1, ++retf){
		//assigning iterators
		decltype(it1) it0, it2 = std::next(it1);
		decltype(retf) retb;
		if (i == 0){
			if (d.is_closed()){
				it0 = ----d.div.end();
				retb = --ret.end();
			} else continue;
		} else{
			it0 = std::prev(it1);
			retb = std::prev(retf);
		}
		//check if we need to snap
		if ((*retf)[0] != -1) continue;
		//calculate back and forward distances to closest corner points
		int iback = (*retb)[retb->size()-2];
		int iforw = (*retf)[1];
		double backdist = (iback>-1) ? d.dist(d.clen[iback], *it1) : 0;
		double forwdist = (iforw>-1) ? d.dist(*it1, d.clen[iforw]) : 0;
		if (backdist > d.snap*d.dist(*it0, *it1)) backdist = 0;
		if (forwdist > d.snap*d.dist(*it1, *it2)) forwdist = 0;
		if (backdist == 0 && forwdist == 0) continue;
		else if (backdist == 0 && forwdist != 0){
			try_forward_move(d, it1, *retb, *retf);
		} else if (backdist != 0 && forwdist == 0){
			try_backward_move(d, it1, *retb, *retf);
		} else if (backdist < forwdist){
			if (!try_backward_move(d, it1, *retb, *retf)){
				try_forward_move(d, it1, *retb, *retf);
			}
		} else {
			if (!try_forward_move(d, it1, *retb, *retf)){
				try_backward_move(d, it1, *retb, *retf);
			}
		}
	}
	return ret;
}

vector<vector<Point*>> final_result(data& d, const vector<vector<int>>& r, HMCont2D::PCollection& pcol){
	vector<vector<Point*>> ret(r.size());
	auto it = d.div.begin();
	for (int i=0; i<r.size(); ++i, ++it){
		Point* p1 = (r[i][0] == -1) ? pcol.add_value(d.get_point(*it)) : d.cont[r[i][0]];
		ret[i].push_back(p1);
		for (int j=1; j<r[i].size()-1; ++j) ret[i].push_back(d.cont[r[i][j]]);
	}
	if (!d.is_closed()) ret.push_back({d.cont.back()});
	return ret;
}
}


vector<vector<Point*>> HMCont2D::Algos::Coarsening(const HMCont2D::Contour& cont, 
		const std::set<Point*>& mandatory,
		HMCont2D::PCollection& pcol,
		double h, double crit_angle, double snap_dist){
	//fill data
	data d;
	d.cont = cont.ordered_points();
	d.clen = {0};
	for (int i=1; i<d.cont.size(); ++i){
		d.clen.push_back(d.clen[i-1] + Point::dist(*d.cont[i], *d.cont[i-1]));
	}
	d.h = h;
	d.hdiamond = tan(crit_angle/2.0);
	d.snap = snap_dist;

	//mandatory points
	for (auto p: mandatory){
		auto ifnd = std::find(d.cont.begin(), d.cont.end(), p) - d.cont.begin();
		if (ifnd < d.cont.size()){
			d.div.push_back(d.clen[ifnd]);
		}
	}
	d.div.sort();
	if (d.div.size()==0 || d.div.front() != 0) d.div.push_front(0);
	if (d.div.back() != d.clen.back()) d.div.push_back(d.clen.back());

	//1. angle points
	add_angle_points(d, crit_angle/2.0);

	//2. initial uniform division
	uniform_division(d);
	
	//3. apply angle restriction
	diamond_rule(d);

	//4. smooth to satisfy len<3*adj_len restriction
	one_third_rule(d);

	//5. snapping
	vector<vector<int>> snres = snapping(d);

	//6. assemble result
	return final_result(d, snres, pcol);
}
