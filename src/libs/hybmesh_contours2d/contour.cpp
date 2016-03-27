#include "contour.hpp"
#include "tree.hpp"
#include "clipper_core.hpp"

using namespace HMCont2D;

Point* Contour::first() const{
	if (size() == 0) return 0;
	else if (size() == 1) return data[0]->pstart;
	else return Edge::PointOrder(*data[0], *data[1])[0];
}

Point* Contour::last() const{
	if (size() == 0) return 0;
	else if (size() == 1) return data[0]->pend;
	else return Edge::PointOrder(*data[size()-2], *data[size()-1])[2];
}

vector<Point*> Contour::ordered_points() const{
	if (size() == 0) return {};
	else if (size() == 1) return {data[0]->pstart, data[0]->pend};
	else {
		vector<Point*> ret {first()};
		for (int i=0; i<size()-1; ++i){
			ret.push_back(Edge::PointOrder(*data[i], *data[i+1])[1]);
		}
		ret.push_back(last());
		return ret;
	}
}

vector<Point*> Contour::unique_points() const{
	auto op = ordered_points();
	if (op.size() == 0) return op;
	vector<Point*> ret {op[0]};
	for (int i=1; i<op.size(); ++i){
		if (*op[i] != *ret.back()) ret.push_back(op[i]);
	}
	return ret;
}

vector<Point*> Contour::corner_points() const{
	auto pnt = unique_points();
	if (pnt.size() == 0) return {};
	bool cl = is_closed();
	vector<Point*> ret;
	if (!cl) ret.push_back(pnt[0]);
	for (int i=0; i<(int)pnt.size()-1; ++i){
		Point *pprev, *p, *pnext;
		if (i == 0){
			if (cl) pprev = pnt[pnt.size() - 2];
			else continue;
		} else pprev = pnt[i-1];
		p = pnt[i];
		pnext = pnt[i+1];
		//build triangle. if its area = 0 -> that is not a corner point
		double ksi = Point::meas_section(*p, *pprev, *pnext);
		if (ksi>=geps*geps){
			ret.push_back(p);
		}
	}
	if (!cl) ret.push_back(pnt.back());
	return ret;
}

vector<Point*> Contour::corner_points1() const{
	auto ret = corner_points();
	if (ret.size()>0 && is_closed()){
		ret.push_back(ret[0]);
	}
	return ret;
}

std::array<Point*, 3> Contour::point_siblings(Point* p) const{
	std::array<Point*, 3> ret {0, 0, 0};
	auto e1 = std::find_if(data.begin(), data.end(), [&p](shared_ptr<Edge> e){ return e->contains(p); });
	if (e1 == data.end()) return ret;
	else { ret[0] = (*e1)->sibling(p); ret[1] = p; }

	auto e2 = std::find_if(e1 + 1, data.end(), [&p](shared_ptr<Edge> e){ return e->contains(p); });
	if (e2 == data.end()){
		assert(!is_closed() && (p==first() || p==last()));
		if (p == first()) std::swap(ret[0], ret[2]);
		return ret;
	} else ret[2] = (*e2)->sibling(p);

	if (is_closed() && (e1==data.begin() && e2 != data.begin()+1)) std::swap(ret[0], ret[2]);

	return ret;
}

std::array<Point*, 3> Contour::point_siblings(int i) const{
	assert((is_closed() && i<size()) || (!is_closed() && i<=size()));
	//get edges
	//e1: next edge
	HMCont2D::Edge *e1 = data[i].get(), *e2;
	//e2: previous edge
	if (i>0 && i<size()){
		e2 = data[i-1].get();
	} else if (is_closed()){
		if (i==0) e2 = data.back().get();
		else e2 = data[i-1].get();
	} else {
		e2 = 0;
		if (i == size()) std::swap(e1,e2);
	}

	//siblings
	if (e1!=0 && e2!=0){
		if (e1->pstart == e2->pstart) return {e2->pend, e1->pstart, e1->pend};
		if (e1->pstart == e2->pend) return {e2->pstart, e1->pstart, e1->pend};
		if (e1->pend == e2->pstart) return {e2->pend, e1->pend, e1->pstart};
		//if e1->pend == e2->pend
		return {e2->pstart, e1->pend, e1->pstart};
	} else if (e1!=0){
		Point* ps = first();
		return {0, ps, e1->sibling(ps)};
	} else {
		Point* ps = last();
		return {e2->sibling(ps), ps, 0};
	}
}

bool Contour::correctly_directed_edge(int i) const{
	if (size() == 1) return true;
	if (is_closed() && size() == 2) return true;
	Point* p2;
	Edge* e2;
	if (i == 0){
		p2 = edge(i)->pend;
		e2 = edge(i+1);
	} else {
		p2 = edge(i)->pstart;
		e2 = edge(i-1);
	}
	return e2->contains(p2);
}

Contour::PInfo Contour::pinfo(Point* p) const{
	auto oi = ordered_info();
	for (auto& i: oi){
		if (i.p == p) return i;
	}
	assert(false);
}
vector<Contour::PInfo> Contour::ordered_info() const{
	vector<PInfo> ret;
	if (size() == 0) return ret;
	assert(is_closed() ? (size() > 2) : true);
	vector<Point*> pts = ordered_points();
	for (int i=0; i<pts.size(); ++i){
		ret.push_back(PInfo());
		auto& a = ret.back();
		a.index = i;
		a.p = pts[i];
		if (i == 0) {
			if (is_closed()){
				a.pprev = pts[pts.size() - 2];
				a.eprev = data.back();
			} else {
				a.pprev = 0;
				a.eprev.reset();
			}
		} else {
			a.pprev = pts[i-1];
			a.eprev = data[i-1];
		}
		if (i == pts.size() - 1){
			if (is_closed()){
				a.pnext = pts[1];
				a.enext = data[0];
			} else {
				a.pnext = 0;
				a.enext.reset();
			}
		} else {
			a.pnext = pts[i+1];
			a.enext = data[i];
		}

	}
	return ret;
}

std::tuple<double, double, int, double, double>
Contour::coord_at(const Point& p) const{
	auto fnd = FindClosestEdge(*this, p);
	int ind = get_index(std::get<0>(fnd));
	auto lens = ELengths(*this);
	double outlen = 0;
	for (int i=0; i<ind; ++i) outlen+=lens[i];

	if (correctly_directed_edge(ind))
		outlen += lens[ind]*std::get<2>(fnd);
	else
		outlen += lens[ind]*(1-std::get<2>(fnd));

	double outw = outlen/std::accumulate(lens.begin(), lens.end(), 0.0);
	return std::make_tuple(outlen, outw, ind, std::get<2>(fnd), std::get<1>(fnd));
}

void Contour::DirectEdges(){
	auto p = ordered_points();
	for (int i=0; i<size(); ++i){
		Point* p0 = p[i];
		Point* p1 = p[i+1];
		if (data[i]->pstart != p0) data[i]->Reverse();
		assert(data[i]->pstart == p0 && data[i]->pend == p1);
	}
}

void Contour::ReallyReverse(){
	if (size() == 1) data[0]->Reverse();
	else{
		Reverse();
		DirectEdges();
	}
}

void Contour::AddLastPoint(Point* p){
	auto np = _generator.allocate();
	np->pstart = last();
	np->pend = p;
	add_value(np);
}

void Contour::RemoveEdge(int i){
	Point *p1 = data[i]->pstart, *p2 = data[i]->pend;
	auto connect = [&](int k){
		if (data[k]->pstart == p2) data[k]->pstart = p1;
		if (data[k]->pend == p2) data[k]->pend = p1;
	};
	if (size()>2){
		if (!is_closed()){
			if (i!=0 && i==size()-1){
				connect(i+1);
				connect(i-1);
			}
		} else {
			connect( (i==0) ? size()-1 : i-1 );
			connect( (i==size()-1) ? 0 : i+1 );
		}
	}
	data.erase(data.begin()+i);
}

void Contour::RemovePoint(const Point* p){
	auto op = ordered_points();
	int edge_to_delete = -1;
	auto change_edge_point = [&](int i, const Point* from, Point* to){
		if (data[i]->pstart == from) data[i]->pstart = to;
		else if (data[i]->pend == from) data[i]->pend = to;
	};
	for (int i=0; i<op.size(); ++i){
		if (op[i] == p){
			edge_to_delete = i;
			if (!is_closed()){
				if (i == op.size() - 1) edge_to_delete = i-1;
				else if (i != 0) change_edge_point(i-1, p, op[i+1]);
			} else {
				if (i == 0) change_edge_point(size()-1, p, op[i+1]);
				else change_edge_point(i-1, p, op[i+1]);
			}
			break;
		}
	}
	if (edge_to_delete>=0){
		data[edge_to_delete]->pstart = 0;
		data[edge_to_delete]->pend = 0;
		data.erase(data.begin()+edge_to_delete);
	}
}

void Contour::RemovePoints(const vector<const Point*>& vp){
	for (auto pnt: vp) RemovePoint(pnt);
}

bool Contour::ForceDirection(bool dir){
	if (Area(*this) < 0){
		if (dir){ Reverse(); return true;}
	} else {
		if (!dir){ Reverse(); return true;}
	}
	return false;
}

bool Contour::IsWithin(const Point& p) const{
	assert(is_closed());
	//use clipper procedure
	Impl::ClipperPath cp(*this);
	return cp.WhereIs(p) == 1;
}

bool Contour::IsWithout(const Point& p) const{
	assert(is_closed());
	//use clipper procedure
	Impl::ClipperPath cp(*this);
	return cp.WhereIs(p) == 0;
}

bool Contour::AllWithout(const vector<Point>& p) const{
	assert(is_closed());
	Impl::ClipperPath cp(*this);
	for (auto& it: p) if (cp.WhereIs(it) != 0) return false;
	return true;
}
int Contour::WhereIs(const Point& p) const{
	assert(is_closed());
	//use clipper procedure
	Impl::ClipperPath cp(*this);
	return cp.WhereIs(p);
}

namespace {
Point inner_point_core(const Contour& c){
	//c is closed inner contour
	//1) build vector from first point of c along the median
	//   of respective angle
	auto tri = c.point_siblings(0);
	//this is a workaround if contour has doubled points
	int iedge1 = 1, iedge2 = c.size()-1;
	if (*tri[2] == *tri[1]){
		for (int i=2; i<=c.size()-2; ++i){
			if (*tri[2] == *tri[1]){
				tri[2] = c.next_point(tri[2]);
				++iedge1;
			} else { break; }
		}
	}
	if (*tri[0] == *tri[1]){
		for (int i=c.size()-2; i>=1; --i){
			if (*tri[0] == *tri[1]){
				tri[0] = c.prev_point(tri[0]);
				--iedge2;
			} else { break; }
		}
	}

	//degenerate cases
	int trieq = 0;
	if (*tri[0] == *tri[1]) trieq+=1;
	if (*tri[0] == *tri[2]) trieq+=2;
	switch (trieq){
		case 1: return Point::Weigh(*tri[0], *tri[2], 0.5);
		case 2: return Point::Weigh(*tri[0], *tri[1], 0.5);
		case 3: return *tri[0];
	}
	
	double area3 = triarea(*tri[0], *tri[1], *tri[2]);
	Vect v;
	if (fabs(area3)<geps*geps){
		//perpendicular to first edge
		v = Vect(-tri[2]->y+tri[1]->y, tri[2]->x-tri[1]->x);
	} else {
		//median point - tri[1]
		Point cp = Point::Weigh(*tri[0], *tri[2], 0.5);
		v = cp - *tri[1];
		if (area3 < 0) v *= -1;
	}
	//2) make this vector longer then the contour size
	BoundingBox box = Contour::BBox(c);
	double len = (box.lenx() + box.leny());
	vecSetLen(v, len);
	//3) find first cross point of vector and contour
	Point a0 = *c.first();
	Point a1 = a0 + v;
	double ksieta[2];
	//using iedge1, iedge2 which are normally equal to 1 and size()-1
	//until first point is doulbled
	for (int i=iedge1; i<iedge2; ++i){
		const Edge* e = c.data[i].get();
		Point b0 = *e->pstart;
		Point b1 = *e->pend;
		SectCrossWRenorm(a0, a1, b0, b1, ksieta);
		if (ksieta[1]>-geps && ksieta[1]<1+geps){
			Point crossp = Point::Weigh(b0, b1, ksieta[1]);
			return Point::Weigh(a0, crossp, 0.5);
		}
	}
	throw std::runtime_error("Failed to find contour inner point");
}
}//inner_point
Point Contour::InnerPoint() const{
	assert(is_closed());
	double a = Area(*this);
	assert(size()>1);
	if (fabs(a)<geps*geps){
		//chose arbitrary point on the edge
		return Point::Weigh(*data[0]->pstart, *data[0]->pend, 0.543);
	} else if (a<0){
		Contour c2;
		c2.data.insert(c2.data.end(), data.rbegin(), data.rend());
		return inner_point_core(c2);
	} else return inner_point_core(*this);
}

std::tuple<bool, Point*>
Contour::GuaranteePoint(const Point& p, PCollection& pcol){
	std::tuple<bool, Point*> ret;
	Point* pc = FindClosestNode(*this, p);
	if ( p == *pc){
		std::get<0>(ret) = false;
		std::get<1>(ret) = pc;
	} else {
		auto ce = FindClosestEdge(*this, p);
		if (ISEQ(std::get<2>(ce), 0)){
			std::get<0>(ret) = false;
			std::get<1>(ret) = std::get<0>(ce)->pstart;
		} else if (ISEQ(std::get<2>(ce), 1)){
			std::get<0>(ret) = false;
			std::get<1>(ret) = std::get<0>(ce)->pend;
		} else {
			Point* p1 = std::get<0>(ce)->pstart;
			Point* p2 = std::get<0>(ce)->pend;
			shared_ptr<Point> pnew(new Point(Point::Weigh(*p1, *p2, std::get<2>(ce))));
			pcol.add_value(pnew);
			int ind = get_index(std::get<0>(ce));
			bool dircorrect = correctly_directed_edge(ind);
			//cannot simply rewrite data to indexed edge because it can be
			//used by another owners.
			RemoveAt({ind});
			auto e1 = std::make_shared<Edge>(p1, pnew.get());
			auto e2 = std::make_shared<Edge>(pnew.get(), p2);
			if (dircorrect) AddAt(ind, {e1, e2});
			else AddAt(ind, {e2, e1});
			std::get<0>(ret) = true;
			std::get<1>(ret) = pnew.get();
		}
	}
	return ret;
}

double Contour::Area(const Contour& c){
	assert(c.is_closed());
	auto order = c.ordered_points();
	double ret = 0;
	auto get_point = [&order](int i){ return order[i]; };
	auto p1 = get_point(0);
	for (int i=1; i<order.size()-2; ++i){
		auto p2 = order[i];
		auto p3 = order[i+1];
		ret += triarea(*p1, *p2, *p3);
	}
	return ret;
};

PCollection Contour::WeightPoints(const Contour& p, vector<double> vw){
	double len = p.length();
	for (auto& v: vw) v*=len;
	return WeightPointsByLen(p, vw);
}

PCollection Contour::WeightPointsByLen(const Contour& p, vector<double> vw){
	PCollection ret;
	auto pseq = p.ordered_points();
	std::sort(vw.begin(), vw.end());
	assert(!ISLOWER(vw[0], 0) && !ISGREATER(vw.back(), p.length()));

	vector<double> elens {0};
	for (auto& e: p.data) elens.push_back(Point::dist(*e->pstart, *e->pend));
	std::partial_sum(elens.begin(), elens.end(), elens.begin());

	int icur=1;
	for (auto w: vw){
		while (icur<pseq.size()-1 && ISEQLOWER(elens[icur], w)) ++icur;
		Point *pprev = pseq[icur-1], *pnext = pseq[icur];
		double t = (w - elens[icur-1])/(elens[icur] - elens[icur-1]);
		ret.add_value( (*pseq[icur-1]) * (1-t) + (*pseq[icur]) * t );
	}
	return ret;
}

vector<double> Contour::EWeights(const Contour& c){
	vector<double> ret {0};
	vector<double> lens = ELengths(c);
	std::copy(lens.begin(), lens.end(), std::back_inserter(ret));
	std::partial_sum(ret.begin(), ret.end(), ret.begin());
	double L = ret.back();
	std::for_each(ret.begin(), ret.end(), [&L](double& x){ x/=L; });
	return ret;
}
