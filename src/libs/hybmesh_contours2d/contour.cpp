#include "contour.hpp"
#include "finder2d.hpp"
#include "treverter2d.hpp"
#include "debug2d.hpp"
using namespace HM2D;
using namespace HM2D::Contour;
namespace hc = HM2D::Contour;

bool Contour::IsContour(const EdgeData& ed){
	for (int i=0; i<ed.size()-1; ++i){
		if (!ed[i]->connected_to(*ed[i+1])) return false;
	}
	return true;
}

shared_ptr<Vertex> Contour::First(const EdgeData& ed){
	if (ed.size() == 0) return 0;
	else if (ed.size() == 1) return ed[0]->first();
	else {
		if (ed[0]->last() == ed[1]->first() ||
		    ed[0]->last() == ed[1]->last())
			return ed[0]->first();
		else return ed[0]->last();
	}
}

shared_ptr<Vertex> Contour::Last(const EdgeData& ed){
	if (ed.size() == 0) return 0;
	else if (ed.size() == 1) return ed[0]->last();
	else {
		auto& e1 = ed.end()[-2];
		auto& e2 = ed.end()[-1];
		if (e2->first() == e1->last() ||
		    e2->first() == e1->first())
			return e2->last();
		else return e2->first();
	}
}

bool hc::IsClosed(const EdgeData& cont){
	if (cont.size() < 2) {
		return false;
	} else if (cont.size() == 2){
		if (cont[0]->first() == cont.back()->first() &&
		    cont[0]->last() == cont.back()->last()) return true;
		if (cont[0]->last() == cont.back()->first() &&
		    cont[0]->first() == cont.back()->last()) return true;
		return false;
	} else return (cont[0]->first() == cont.back()->last() ||
	               cont[0]->last() == cont.back()->first() ||
	               cont[0]->first() == cont.back()->first() ||
	               cont[0]->last() == cont.back()->last());
}
bool hc::IsOpen(const EdgeData& cont){
	return !IsClosed(cont);
}

VertexData hc::OrderedPoints(const EdgeData& cont){
	if (cont.size() == 0) return {};
	else if (cont.size() == 1) return {cont[0]->first(), cont[0]->last()};
	else {
		VertexData ret {First(cont)};
		for (int i=0; i<cont.size()-1; ++i){
			auto& e1 = cont[i];
			if (e1->first() != ret.back()) ret.push_back(e1->first());
			else ret.push_back(e1->last());
		}
		ret.push_back(Last(cont));
		return ret;
	}
}
VertexData hc::OrderedPoints1(const EdgeData& cont){
	auto ret = OrderedPoints(cont);
	if (ret[0] == ret.back()) ret.resize(ret.size()-1);
	return ret;
}

VertexData hc::UniquePoints(const EdgeData& ed){
	auto op = OrderedPoints(ed);
	if (op.size() == 0) return op;
	VertexData ret {op[0]};
	for (int i=1; i<op.size(); ++i){
		if (*op[i] != *ret.back()) ret.push_back(op[i]);
	}
	return ret;
}
VertexData hc::UniquePoints1(const EdgeData& ed){
	auto ret = UniquePoints(ed);
	if (ret[0] == ret.back()) ret.pop_back();
	return ret;
}

VertexData hc::CornerPoints(const EdgeData& ed){
	auto pnt = UniquePoints(ed);
	if (pnt.size() == 0) return {};
	VertexData ret;
	bool cl = IsClosed(ed);
	if (!cl) ret.push_back(pnt[0]);
	for (int i=0; i<(int)pnt.size()-1; ++i){
		Point *pprev, *p, *pnext;
		if (i == 0){
			if (cl) pprev = pnt[pnt.size() - 2].get();
			else continue;
		} else pprev = pnt[i-1].get();
		p = pnt[i].get();
		pnext = pnt[i+1].get();
		double ksi = Point::meas_section(*p, *pprev, *pnext);
		if (ksi>=geps*geps){
			ret.push_back(pnt[i]);
		}
	}
	if (!cl) ret.push_back(pnt.back());
	return ret;
}

VertexData Contour::CornerPoints1(const EdgeData& ed){
	auto ret = CornerPoints(ed);
	if (ret.size()>0 && Contour::IsClosed(ed)){
		ret.push_back(ret[0]);
	}
	return ret;
}

bool Contour::CorrectlyDirectedEdge(const EdgeData& cont, int i){
	if (cont.size() == 1) return true;
	if (IsClosed(cont) && cont.size() == 2) return true;
	Point* p2;
	Edge* e2;
	if (i == 0){
		p2 = cont[i]->last().get();
		e2 = cont[i+1].get();
	} else {
		p2 = cont[i]->first().get();
		e2 = cont[i-1].get();
	}
	return e2->first().get() == p2 || e2->last().get() == p2;
}

PInfoR Contour::PInfo(const EdgeData& ed, const Point* p){
	auto oi = OrderedInfo(ed);
	for (auto& i: oi){
		if (i.p.get() == p) return i;
	}
	assert(false);
}
vector<PInfoR> Contour::OrderedInfo(const EdgeData& ed){
	vector<PInfoR> ret;
	if (ed.size() == 0) return ret;
	bool is_closed = IsClosed(ed);
	assert(is_closed ? (ed.size() > 2) : true);
	VertexData pts = OrderedPoints(ed);
	for (int i=0; i<pts.size(); ++i){
		ret.push_back(PInfoR());
		auto& a = ret.back();
		a.index = i;
		a.p = pts[i];
		if (i == 0) {
			if (is_closed){
				a.pprev = pts[pts.size() - 2];
				a.eprev = ed.back();
			} else {
				a.pprev = 0;
				a.eprev.reset();
			}
		} else {
			a.pprev = pts[i-1];
			a.eprev = ed[i-1];
		}
		if (i == pts.size() - 1){
			if (is_closed){
				a.pnext = pts[1];
				a.enext = ed[0];
			} else {
				a.pnext = 0;
				a.enext.reset();
			}
		} else {
			a.pnext = pts[i+1];
			a.enext = ed[i];
		}

	}
	return ret;
}

std::tuple<double, double, int, double, double>
hc::CoordAt(const EdgeData& cont, const Point& p){
	auto fnd = HM2D::Finder::ClosestEdge(cont, p);
	int ind = std::get<0>(fnd);
	assert(ind >= 0);
	auto lens = ELengths(cont);
	double outlen = std::accumulate(lens.begin(), lens.begin()+ind, 0.0);

	if (CorrectlyDirectedEdge(cont, ind))
		outlen += lens[ind]*std::get<2>(fnd);
	else
		outlen += lens[ind]*(1-std::get<2>(fnd));

	double outw = outlen/std::accumulate(lens.begin(), lens.end(), 0.0);
	return std::make_tuple(outlen, outw, ind, std::get<2>(fnd), std::get<1>(fnd));
}

void hc::Reverse(EdgeData& ed){ std::reverse(ed.begin(), ed.end()); }

void Contour::AddLastPoint(EdgeData& to, std::shared_ptr<Vertex> p){
	auto p0 = Last(to);
	to.push_back(std::make_shared<Edge>(p0, p));
}

namespace {
Point midp(const VertexData& up){
	Point ret(*up[0]);
	for (int i=1; i<up.size(); ++i){ ret += *up[i]; }
	return ret/(up.size());
}
bool inside_tri(Point p, Point t1, Point t2, Point t3){
	p-=t1; t2-=t1; t3-=t1;
	double modj = (t3.y*t2.x - t3.x*t2.y);
	assert(modj>geps*geps);
	double ksi = ( t3.y*(p.x) - t3.x*(p.y));
	double eta = (-t2.y*(p.x) + t2.x*(p.y));
	return (ksi>0 && ksi<modj && eta>0 && eta<modj-ksi);
}
Point inner_point_core(const EdgeData& c){
	auto up = UniquePoints1(c);
	//triangle
	if (up.size() == 3){ return midp(up); }
	//find best convex triangle
	vector<double> triareas(up.size());
	for (int i=0; i<up.size(); ++i){
		int ip = (i==0) ? up.size()-1 : i-1;
		int in = (i==up.size()-1) ? 0 : i+1;
		triareas[i] = triarea(*up[ip], *up[i], *up[in]);
	}
	//if convex polygon return average point
	if (std::all_of(triareas.begin(), triareas.end(), [](double a){ return a>-geps*geps; })){
		return midp(up);
	}
	//else find closest point to best convex vertex lying within its triangle
	Point pbest;
	double mbest=std::numeric_limits<double>::max();
	int iconv = std::max_element(triareas.begin(), triareas.end())-triareas.begin();
	if (iconv == 0) std::rotate(up.begin(), up.end()-1, up.end());
	else std::rotate(up.begin(), up.begin()+iconv-1, up.end());
	for (int i=3; i<up.size(); ++i){
		if (inside_tri(*up[i], *up[0], *up[1], *up[2])){
			double m = Point::meas(*up[1], *up[i]);
			if (m<mbest){ mbest = m; pbest = *up[i]; }
		}
	}

	if (mbest == std::numeric_limits<double>::max()) return Point::Weigh(*up[0], *up[2], 0.5);
	else return Point::Weigh(pbest, *up[1], 0.5);
}
}//inner_point

Point Contour::InnerPoint(const EdgeData& ed){
	assert(Contour::IsClosed(ed));
	if (ed.size() == 2) return ed[0]->center();
	assert(ed.size()>2);
	Contour::R::Clockwise rr(ed, false);
	return inner_point_core(ed);
}

std::tuple<bool, shared_ptr<Vertex>>
Contour::GuaranteePoint(EdgeData& ed, const Point& p){
	std::tuple<bool, shared_ptr<Vertex>> ret;

	auto ce = HM2D::Finder::ClosestEdge(ed, p);
	if (std::get<0>(ce)<0){
		std::get<0>(ret) = false;
		return ret;
	}
	Edge* e = ed[std::get<0>(ce)].get();
	double elen = e->length();
	double len2 = std::get<2>(ce)*elen;
	if (ISZERO(len2)){
		std::get<0>(ret) = false;
		std::get<1>(ret) = e->first();
	} else if (ISZERO(elen-len2)){
		std::get<0>(ret) = false;
		std::get<1>(ret) = e->last();
	} else {
		auto pnew = std::make_shared<Vertex>(Point::Weigh(
			*e->first(), *e->last(), std::get<2>(ce)));
		auto e1 = std::make_shared<Edge>(*e);
		auto e2 = std::make_shared<Edge>(*e);
		e1->vertices[1] = pnew;
		e2->vertices[0] = pnew;
		if (CorrectlyDirectedEdge(ed, std::get<0>(ce))){
			ed[std::get<0>(ce)] = e1;
			ed.insert(ed.begin()+std::get<0>(ce)+1, e2);
		} else {
			ed[std::get<0>(ce)] = e2;
			ed.insert(ed.begin()+std::get<0>(ce)+1, e1);
		}
		std::get<0>(ret) = true;
		std::get<1>(ret) = pnew;
	}
	return ret;
}

double hc::Area(const EdgeData& c){
	if (c.size() == 0) return 0;
	assert(IsClosed(c));
	auto order = OrderedPoints(c);
	double ret = 0;
	auto p1 = order[0];
	for (int i=1; i<order.size()-2; ++i){
		auto p2 = order[i];
		auto p3 = order[i+1];
		ret += triarea(*p1, *p2, *p3);
	}
	return ret;
};

Point hc::WeightPoint(const EdgeData& ed, double w){
	return WeightPoints(ed, vector<double>{w})[0];
}

vector<Point> hc::WeightPoints(const EdgeData& p, vector<double> vw){
	double len = Length(p);
	for (auto& v: vw) v*=len;
	return WeightPointsByLen(p, vw);
}

vector<Point> Contour::WeightPointsByLen(const EdgeData& p, vector<double> vw){
	vector<Point> ret;
	auto pseq = OrderedPoints(p);
	std::sort(vw.begin(), vw.end());
	assert(!ISLOWER(vw[0], 0) && !ISGREATER(vw.back(), Length(p)));
	
	vector<double> elens {0};
	for (auto& e: p) elens.push_back(Point::dist(*e->first(), *e->last()));
	std::partial_sum(elens.begin(), elens.end(), elens.begin());

	int icur=1;
	for (auto w: vw){
		while (icur<pseq.size()-1 && ISEQLOWER(elens[icur], w)) ++icur;
		Point *pprev = pseq[icur-1].get(), *pnext = pseq[icur].get();
		double t = (w - elens[icur-1])/(elens[icur] - elens[icur-1]);
		ret.push_back( (*pseq[icur-1]) * (1-t) + (*pseq[icur]) * t );
	}
	return ret;
}

vector<double> Contour::EWeights(const EdgeData& c){
	vector<double> ret {0};
	vector<double> lens = ELengths(c);
	std::copy(lens.begin(), lens.end(), std::back_inserter(ret));
	std::partial_sum(ret.begin(), ret.end(), ret.begin());
	double L = ret.back();
	std::for_each(ret.begin(), ret.end(), [&L](double& x){ x/=L; });
	return ret;
}

void hc::Connect(EdgeData& to, const EdgeData& from){
	auto self0 = First(to), self1 = Last(to);
	auto target0 = First(from), target1 = Last(from);
	//choosing option for unition
	if (to.size() == 0 || from.size() == 0 ) goto COPY12;
	else if (self0 == self1 || target0 == target1) goto THROW;
	else if (from.size() == 1 &&
		(from[0]->first() == self1 || from[0]->last() == self1)) goto COPY12;
	else if (from.size() == 1 &&
		(from[0]->first() == self0 || from[0]->last() == self0)) goto COPY03;
	//try to add new contour to the end of current
	else if (self1 == target0) goto COPY12;
	else if (self1 == target1) goto NEED_REVERSE;
	//if failed try to add before the start
	else if (self0 == target1) goto COPY03;
	else if (self0 == target0) goto NEED_REVERSE;
	else goto THROW;

COPY03:
	{
		to.insert(to.begin(), from.begin(), from.end());
		return;
	}
COPY12:
	{
		to.insert(to.end(), from.begin(), from.end());
		return;
	}
NEED_REVERSE:
	{
		EdgeData tmp(from);
		Reverse(tmp);
		return Connect(to, tmp);
	}
THROW:
	{
		throw std::runtime_error("Impossible to unite non-connected contours");
	}
}

void hc::SplitEdge(EdgeData& cont, int iedge, const vector<Point>& pts){
	if (pts.size() == 0) return;
	assert(IsContour(cont));
	bool iscor = CorrectlyDirectedEdge(cont, iedge);
	if (!iscor) cont[iedge]->reverse();

	VertexData pp(pts.size()+2);
	pp[0] = cont[iedge]->vertices[0];
	pp.back() = cont[iedge]->vertices[1];
	for (int i=1; i<pp.size()-1; ++i){
		pp[i] = std::make_shared<Vertex>(pts[i-1]);
	}

	cont[iedge]->vertices[1] = pp[1];
	EdgeData newed(pts.size());
	for (int i=0; i<pts.size(); ++i){
		newed[i] = std::make_shared<Edge>(*cont[iedge]);
		newed[i]->vertices[0] = pp[i+1];
		newed[i]->vertices[1] = pp[i+2];
	}
	cont.insert(cont.begin()+iedge+1, newed.begin(), newed.end());

	if (!iscor){
		for (auto it=cont.begin()+iedge; it!=cont.begin()+iedge+pts.size()+1; ++it){
			(*it)->reverse();
		}
	}
};
