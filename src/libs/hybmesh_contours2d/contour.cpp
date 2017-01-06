#include "contour.hpp"
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
	if (cont.size() < 2) return false;
	return (cont[0]->first() == cont.back()->last() ||
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

VertexData hc::UniquePoints(const EdgeData& ed){
	auto op = OrderedPoints(ed);
	if (op.size() == 0) return op;
	VertexData ret {op[0]};
	for (int i=1; i<op.size(); ++i){
		if (*op[i] != *ret.back()) ret.push_back(op[i]);
	}
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
		e2 = cont[i].get();
	} else {
		p2 = cont[i]->first().get();
		e2 = cont[i-1].get();
	}
	return e2->first().get() == p2 || e2->first().get() == p2;
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
	auto fnd = FindClosestEdge(cont, p);
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

int Contour::WhereIs(const EdgeData& ed, const Point& p){
	//whereis for contour trees uses this routine
	//passing all tree edges as 'ed'.
	//Hence 'ed' is not always a closed contour.
	//assert(IsClosed(ed));
	auto bbox = BBox(ed, 0);
	if (bbox.whereis(p) == OUTSIDE) return OUTSIDE;

	//calculate number of crosses between ed and [p, pout]
	//where pout is random outside point
	Point pout(bbox.xmax + 1.56749, bbox.ymax + 1.06574);

	for (int tries=0; tries<100; ++tries){
		double ksieta[2];
		int ncrosses = 0;
		for (auto& e: ed){
			SectCross(p, pout, *e->first(), *e->last(), ksieta);
			if (ISIN_NN(ksieta[0], 0, 1) && ISIN_NN(ksieta[1], 0, 1)){
				++ncrosses;
			} else if (ISEQ(ksieta[0], 0)){
			       return BOUND;
			} else if (ISEQ(ksieta[1], 0) || ISEQ(ksieta[1], 1)){
				//[p, pout] crosses one of ed vertices.
				//This is ambiguous so we change pout and try again
				goto NEXTTRY;
			}
		}
		if (ncrosses % 2 == 0) return OUTSIDE;
		else return INSIDE;
NEXTTRY:
		pout+=Point(-0.4483, 0.0342);
	}

	throw std::runtime_error("failed to detect point-contour relation");
}

namespace {
Point inner_point_core(const EdgeData& c){
	auto up = UniquePoints(c);
	assert(up.size()>=3);
	//c is closed inner contour
	//1) build vector from first point of c along the median
	//   of respective angle
	std::array<const Point*, 3> tri {up.back().get(), up[0].get(), up[1].get()};
	//this is a workaround if contour has doubled points
	int iedge1 = 1, iedge2 = up.size()-1;

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
	BoundingBox box = HM2D::BBox(up, 0);
	double len = (box.lenx() + box.leny());
	vecSetLen(v, len);
	//3) find first cross point of vector and contour
	Point a0 = *up[0];
	Point a1 = a0 + v;
	double ksieta[2];
	//using iedge1, iedge2 which are normally equal to 1 and size()-1
	//until doubled points are found
	Point best_cross;
	double minksi = std::numeric_limits<double>::max();
	int ilimit=iedge1 - 1;
	//looking for cross going forward
	for (int i=iedge1; i<iedge2; ++i){
		Point b0 = *up[i];
		Point b1 = *up[i+1];
		SectCrossWRenorm(a0, a1, b0, b1, ksieta);
		if (ksieta[0]>geps && ksieta[1]>-geps && ksieta[1]<1+geps){
			minksi = ksieta[0];
			best_cross = Point::Weigh(b0, b1, ksieta[1]);
			ilimit = i;
			break;
		}
	}
	//looking for cross going backward
	for (int i=iedge2-1; i>ilimit; --i){
		Point b0 = *up[i];
		Point b1 = *up[i+1];
		SectCrossWRenorm(a0, a1, b0, b1, ksieta);
		if (ksieta[0]>geps && ksieta[1]>-geps && ksieta[1]<1+geps){
			if (ksieta[0] < minksi){
				minksi = ksieta[0];
				best_cross = Point::Weigh(b0, b1, ksieta[1]);
			}
			break;
		}
	}
	
	assert(minksi < std::numeric_limits<double>::max() - 1.0);
	return Point::Weigh(a0, best_cross, 0.5);
}
}//inner_point

Point Contour::InnerPoint(const EdgeData& ed){
	assert(Contour::IsClosed(ed));
	assert(ed.size()>1);
	double a = Contour::Area(ed);
	if (fabs(a)<geps*geps){
		//chose arbitrary point on the edge
		return Point::Weigh(*ed[0]->first(), *ed[0]->last(), 0.543);
	} else if (a<0){
		EdgeData c2 = ed;
		Reverse(c2);
		return inner_point_core(c2);
	} else return inner_point_core(ed);
}

std::tuple<bool, shared_ptr<Vertex>>
Contour::GuaranteePoint(EdgeData& ed, const Point& p){
	std::tuple<bool, shared_ptr<Vertex>> ret;

	auto ce = FindClosestEdge(ed, p);
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

	COPY03:{
		to.insert(to.begin(), from.begin(), from.end());
		return;
	}
	COPY12:{
		to.insert(to.end(), from.begin(), from.end());
		return;
	}
	NEED_REVERSE:{
		EdgeData tmp(from);
		Reverse(tmp);
		return Connect(to, tmp);
	}
	THROW:{
		throw std::runtime_error("Impossible to unite non-connected contours");
	}
}
