#include "contour.hpp"
#include "finder2d.hpp"
#include "treverter2d.hpp"
#include "debug2d.hpp"
using namespace HM2D;
using namespace HM2D::Contour;

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

bool Contour::IsClosed(const EdgeData& cont){
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
bool Contour::IsOpen(const EdgeData& cont){
	return !IsClosed(cont);
}

VertexData Contour::OrderedPoints(const EdgeData& cont){
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
VertexData Contour::UniquePoints(const EdgeData& ed){
	auto op = OrderedPoints(ed);
	if (op.size() == 0) return op;
	VertexData ret {op[0]};
	for (int i=1; i<op.size(); ++i){
		if (*op[i] != *ret.back()) ret.push_back(op[i]);
		else ret.back() = op[i];
	}
	return ret;
}
namespace{
bool is_corner(const Point& p, const Point& pprev, const Point& pnext){
	double ksi = Point::meas_section(p, pprev, pnext);
	return ksi>=geps*geps;
}
}
VertexData Contour::CornerPoints(const EdgeData& ed){
	auto pnt = UniquePoints(ed);
	if (pnt.size() == 0) return {};
	VertexData ret;
	bool cl = IsClosed(ed);
	if (!cl) ret.push_back(pnt[0]);
	for (int i=0; i<(int)pnt.size()-1; ++i){
		Point *pprev, *p, *pnext;
		if (!cl && i==0) continue;

		pprev = (i > 0) ? pnt[i-1].get() : pnt[pnt.size()-2].get();
		p = pnt[i].get();
		pnext = pnt[i+1].get();

		if (is_corner(*p, *pprev, *pnext)) ret.push_back(pnt[i]);
	}
	if (!cl) ret.push_back(pnt.back());
	else if (ret.size() > 0) ret.push_back(ret[0]);
	return ret;
}

namespace{
EdgeData nozeros(const EdgeData& ed){
	EdgeData ret;
	double sumlen=0;
	for (auto i=0; i<ed.size(); ++i){
		sumlen += ed[i]->length();
		if (sumlen > geps){
			sumlen = 0;
			ret.push_back(ed[i]);
		}
	}
	return ret;
}
VertexData signpoints_open(const EdgeData& ed){
	//remove zero sized edges
	EdgeData nz = nozeros(ed);
	if (nz.size() == 0) return {};
	//add first
	VertexData ret {nz[0]->first()};
	//add middle points
	for (int i=1; i<nz.size(); ++i){
		//by boundary type
		if (nz[i-1]->boundary_type != nz[i]->boundary_type){
			ret.push_back(nz[i]->first());
			continue;
		}
		//by angle
		auto pprev = nz[i-1]->pfirst();
		auto pcur = nz[i]->pfirst();
		auto pnext = (i==nz.size()-1) ? nz.back()->plast() : nz[i+1]->pfirst();
		if (is_corner(*pcur, *pprev, *pnext)) ret.push_back(nz[i]->first());
	}
	//add last
	ret.push_back(nz.back()->last());
	return ret;
}

VertexData signpoints_closed(const EdgeData& ed){
	//remove zero sized edges
	EdgeData nz = nozeros(ed);
	if (nz.size() < 2) return {};
	if (nz.size() == 2) return {nz[0]->first(), nz[0]->last(), nz[0]->first()};
	VertexData ret;
	//add middle points
	for (int i=0; i<nz.size(); ++i){
		auto eprev = (i==0) ? nz.back() : nz[i-1];
		auto enext = (i==nz.size()-1) ? nz[0] : nz[i+1];
		//by boundary type
		if (eprev->boundary_type != nz[i]->boundary_type){
			ret.push_back(nz[i]->first());
			continue;
		}
		//by angle
		auto pprev = eprev->pfirst();
		auto pcur = nz[i]->pfirst();
		auto pnext = enext->pfirst();
		if (is_corner(*pcur, *pprev, *pnext)) ret.push_back(nz[i]->first());
	}
	//add last
	ret.push_back(ret[0]);
	return ret;
}
}

VertexData Contour::SignificantPoints(const EdgeData& ed){
	Contour::R::ReallyDirect rd(ed);
	if (IsClosed(ed)) return signpoints_closed(ed);
	else return signpoints_open(ed);
}


VertexData Contour::OrderedPoints1(const EdgeData& cont){
	auto ret = OrderedPoints(cont);
	if (ret.size()>0 && ret[0] == ret.back()){
		ret.resize(ret.size()-1);
	}
	return ret;
}
VertexData Contour::CornerPoints1(const EdgeData& ed){
	auto ret = CornerPoints(ed);
	if (ret.size()>0 && ret[0] == ret.back()){
		ret.resize(ret.size()-1);
	}
	return ret;
}
VertexData Contour::SignificantPoints1(const EdgeData& ed){
	auto ret = SignificantPoints(ed);
	if (ret.size()>0 && ret[0] == ret.back()){
		ret.resize(ret.size()-1);
	}
	return ret;
}
VertexData Contour::UniquePoints1(const EdgeData& ed){
	auto ret = UniquePoints(ed);
	if (ret.size()>0 && ret[0] == ret.back()){
		ret.resize(ret.size()-1);
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
Contour::CoordAt(const EdgeData& cont, const Point& p){
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
	//else find closest vertex lying within concave triangle
	Point pbest;
	double mbest=std::numeric_limits<double>::max();
	int iconv = std::max_element(triareas.begin(), triareas.end())-triareas.begin();
	if (iconv == 0) std::rotate(up.begin(), up.end()-1, up.end());
	else std::rotate(up.begin(), up.begin()+iconv-1, up.end());
	//using lesser triangle to exclude cases with contacts to 0-2 lines
	Point nd0 = Point::Weigh(*up[0], *up[1], 0.1);
	Point nd2 = Point::Weigh(*up[1], *up[2], 0.9);
	for (int i=3; i<up.size(); ++i){
		if (inside_tri(*up[i], nd0, *up[1], nd2)){
			double m = Point::meas(*up[1], *up[i]);
			if (m<mbest){ mbest = m; pbest = *up[i]; }
		}
	}

	if (mbest == std::numeric_limits<double>::max()) return Point::Weigh(nd0, nd2, 0.5);
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

double Contour::Area(const EdgeData& c){
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

Point Contour::WeightPoint(const EdgeData& ed, double w){
	return WeightPoints(ed, vector<double>{w})[0];
}

vector<Point> Contour::WeightPoints(const EdgeData& p, vector<double> vw){
	double len = Length(p);
	for (auto& v: vw) v*=len;
	return WeightPointsByLen(p, vw);
}

vector<Point> Contour::WeightPointsByLen(const EdgeData& p, vector<double> vw){
	vector<Point> ret;
	if (vw.size() == 0) return ret;
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

double Contour::Length(const EdgeData& ed){
	double ret = 0;
	for (auto& e: ed) ret += e->length();
	return ret;
}
vector<double> Contour::ELengths(const EdgeData& ed){
	vector<double> ret;
	for (auto& e: ed) ret.push_back(e->length());
	return ret;
}

vector<int> Contour::BTypesFromWeights(const EdgeData& cont, const vector<double>& w){
	assert(Contour::IsContour(cont));
	vector<int> ret(w.size());
	vector<double> eweights = HM2D::Contour::EWeights(cont);

	auto ebegin = eweights.begin();
	for (int i=0; i<w.size(); ++i){
		auto fnd = std::upper_bound(ebegin, eweights.end(), w[i]);
		if (fnd != eweights.begin()) --fnd;
		ebegin = fnd;
		ret[i] = cont[fnd - eweights.begin()]->boundary_type;
	}

	return ret;
}
