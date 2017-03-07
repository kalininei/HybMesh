#include "buildcont.hpp"
#include "spmat.hpp"
#include "partition01.hpp"
#include "contour.hpp"
#include "treverter2d.hpp"
#include "assemble2d.hpp"
#include "modcont.hpp"
#include "debug2d.hpp"

using namespace HM2D;
using namespace HM2D::Contour;
namespace cns = HM2D::Contour::Constructor;
namespace ens = HM2D::ECol::Constructor;

EdgeData cns::Circle(int N, double rad, Point cnt){
	Point poc = cnt + Point(rad, 0);
	return Circle(N, cnt, poc);
};

EdgeData cns::Circle(int N, Point cnt, Point poc){
	EdgeData res;
	Point p = poc-cnt;
	double rad = Point::dist(poc, cnt);
	double angle0 = atan2(p.y, p.x);
	VertexData pdata;
	for (int i=0; i<N; ++i){
		double angle = 2*M_PI/N*i + angle0;
		pdata.emplace_back(new Vertex(rad*cos(angle), rad*sin(angle)));
		*pdata.back() += cnt;
		if (i>0) res.emplace_back(new Edge(pdata[i-1], pdata[i]));
	}
	res.emplace_back(new Edge(pdata[N-1], pdata[0]));
	return res;
}

EdgeData cns::FromPoints(const vector<double>& pnt, bool force_closed){
	vector<Point> p;
	for (int i=0; i<pnt.size(); i+=2) p.emplace_back(pnt[i], pnt[i+1]);
	return FromPoints(p, force_closed);
}

EdgeData cns::FromPoints(const vector<Point>& pnt, bool force_closed){
	if (pnt.size() == 0) return HM2D::EdgeData();
	VertexData p;
	for (int i=0; i<(int)pnt.size(); ++i) p.emplace_back(new Vertex(pnt[i]));
	if (*p[0] == *p.back()) p.back() = p[0];
	return Contour::Assembler::Contour1(p, force_closed);
}

EdgeData cns::CutContour(const EdgeData& cont, const Point& pstart, int direction, double len){
	if (direction == -1){
		R::ReallyRevert rr(cont);
		return CutContour(cont, pstart, 1, len);
	}
	EdgeData other;
	DeepCopy(cont, other);
	R::ReallyDirect::Permanent(other);
	auto ps = std::get<1>(Algos::GuaranteePoint(other, pstart));
	if (IsClosed(other)){
		R::ForceFirst::Permanent(other, pstart);
	} else {
		int ind = std::find_if(other.begin(), other.end(),
			[&ps](shared_ptr<Edge> e){ return e->first() == ps; })
			- other.begin();
		assert(ind < other.size());
		other.erase(other.begin(), other.begin() + ind);
	}
	assert(len<=Length(other));
	Point p2 = Contour::WeightPointsByLen(other, {len})[0];
	auto ps2 = std::get<1>(Algos::GuaranteePoint(other, p2));
	return Assembler::ShrinkContour(other, ps.get(), ps2.get());
}

EdgeData cns::CutContour(const EdgeData& cont, const Point& pstart, const Point& pend){
	EdgeData ret;
	DeepCopy(cont, ret);
	auto p1 = std::get<1>(Algos::GuaranteePoint(ret, pstart));
	auto p2 = std::get<1>(Algos::GuaranteePoint(ret, pend));
	ret = Assembler::ShrinkContour(ret, p1.get(), p2.get());
	if (First(ret) != p1) R::ReallyRevert::Permanent(ret);
	return ret;
}

EdgeData ens::FromRaw(int npnt, int neds, double* pnt, int* eds){
	VertexData pdata;
	for (int i=0; i<npnt; ++i) pdata.push_back(std::make_shared<Vertex>(pnt[2*i], pnt[2*i+1]));
	EdgeData ret;
	for (int i=0; i<neds; ++i){
		ret.push_back(std::make_shared<Edge>(
			pdata[eds[2*i]], pdata[eds[2*i+1]]));
	}
	return ret;
}

namespace {
//A*k0 + B*k1 + C
std::array<double, 3> second_der(double x0, double y0, double x1, double y1, double t){
	double mult = 2.0/(x1 - x0)/(x1 - x0);
	return {(3*t-2)*(x1-x0)*mult, (3*t-1)*(x1 - x0)*mult, (3-6*t)*(y1-y0)*mult};
}

//A*k0 + B*k1 + C*k2 + D
std::array<double, 4> second_der(double x0, double y0, double x1, double y1, double x2, double y2){
	auto p1 = second_der(x0, y0, x1, y1, 1.0);
	auto p2 = second_der(x1, y1, x2, y2, 0.0);
	return {p1[0], p1[1] - p2[0], -p2[1], p1[2]-p2[2]};
}

HMMath::Mat spline_k_internal(const vector<double>& x, const vector<double>& y, vector<double>& rhs){
	HMMath::Mat mat;
	for (int i=1; i<x.size()-1; ++i){
		std::array<double, 4> p = second_der(x[i-1], y[i-1], x[i], y[i], x[i+1], y[i+1]);
		mat.set(i, i-1, p[0]);
		mat.set(i, i, p[1]);
		mat.set(i, i+1, p[2]);
		rhs[i] = -p[3];
	}
	return mat;
}

vector<double> spline_k_natural(const vector<double>& x, const vector<double>& y){
	vector<double> rhs(x.size(), 0);
	auto mat = spline_k_internal(x, y, rhs);
	//boundary conditions
	std::array<double, 3> p1 = second_der(x[0], y[0], x[1], y[1], 0.0);
	mat.set(0, 0, p1[0]);
	mat.set(0, 1, p1[1]);
	rhs[0] = -p1[2];
	std::array<double, 3> p2 = second_der(x[x.size()-2], y[y.size()-2], x.back(), y.back(), 1.0);
	mat.set(x.size()-1, x.size()-2, p2[0]);
	mat.set(x.size()-1, x.size()-1, p2[1]);
	rhs.back() = -p2[2];
	//solution
	auto slv = HMMath::MatSolve::Factory(mat);
	vector<double> ans(x.size(), 0);
	slv->Solve(rhs, ans);
	return ans;
}

vector<double> spline_k_periodic(const vector<double>& x, const vector<double>& y){
	vector<double> rhs(x.size(), 0);
	auto mat = spline_k_internal(x, y, rhs);
	//boundary conditions
	std::array<double, 4> p = second_der(x[0] - (x.back() - x[x.size()-2]), y.back(), x[0], y[0], x[1], y[1]);
	mat.set(0, x.size()-2, p[0]);
	mat.set(0, 0, p[1]);
	mat.set(0, 1, p[2]);
	rhs[0] = - p[3];
	mat.set(x.size()-1, x.size()-1, 1.0);
	mat.set(x.size()-1, 0, -1.0);
	//solution
	auto slv = HMMath::MatSolve::Factory(mat);
	vector<double> ans(x.size(), 0);
	slv->Solve(rhs, ans);
	return ans;
}

std::array<double, 4> to_spline_segment(double x0, double x1, double y0, double y1, double k1, double k2){
	double a = k1*(x1-x0) - (y1-y0);
	double b = -k2*(x1-x0) + (y1-y0);
	return {a-b, -2*a+b, -y0+y1+a, y0};
}

};

EdgeData cns::Spline(const vector<Point>& pnt1, int nedges, bool force_closed){
	//force closed option
	if (force_closed && pnt1[0] != pnt1.back()){
		auto pnew = pnt1;
		pnew.push_back(pnt1[0]);
		return Spline(pnew, nedges, false);
	}
	vector<const Point*> pnt;
	for (auto& p: pnt1) pnt.push_back(&p);
	if (*pnt[0] == *pnt.back()) pnt.back() = pnt[0];

	int nseg = pnt.size() - 1;
	vector<double> double_nsum;
	vector<double> vlen(1, 0.0);
	for (int i=0; i<nseg; ++i){
		double len = Point::dist(*pnt[i+1], *pnt[i]);
		vlen.push_back(vlen.back() + len);
		double_nsum.push_back(len);
	}
	for (int i=0; i<nseg; ++i) double_nsum[i] *= (nedges/vlen.back());
	vector<int> nsum = HMMath::RoundVector(double_nsum, vector<int>(nseg, 1));

	std::vector<double> xvec, yvec;
	for (auto p: pnt) {xvec.push_back(p->x); yvec.push_back(p->y); }
	std::vector<double> kx = (pnt[0] != pnt.back()) ? spline_k_natural(vlen, xvec)
		                                        : spline_k_periodic(vlen, xvec);
	std::vector<double> ky = (pnt[0] != pnt.back()) ? spline_k_natural(vlen, yvec)
		                                        : spline_k_periodic(vlen, yvec);
	std::vector<std::array<double, 4>> xcubic(nseg), ycubic(nseg);
	for (int i=0; i<nseg; ++i){
		xcubic[i] = to_spline_segment(vlen[i], vlen[i+1], xvec[i], xvec[i+1], kx[i], kx[i+1]);
		ycubic[i] = to_spline_segment(vlen[i], vlen[i+1], yvec[i], yvec[i+1], ky[i], ky[i+1]);
	}
	vector<Point> pp;
	for (int i=0; i<nseg; ++i){
		pp.push_back(*pnt[i]);
		auto& fx = xcubic[i];
		auto& fy = ycubic[i];
		for (int j=1; j<nsum[i]; ++j){
			double t = (double)j/(nsum[i]);
			double x = fx[0]*t*t*t + fx[1]*t*t + fx[2]*t + fx[3];
			double y = fy[0]*t*t*t + fy[1]*t*t + fy[2]*t + fy[3];
			pp.push_back(Point(x, y));
		}
	}
	pp.push_back(*pnt.back());

	return cns::FromPoints(pp, false);
}

EdgeData cns::FromContours(const vector<EdgeData>& conts,
	bool last_close, bool fullshift, std::set<int> fixed){
	//checks
	for (int i=0; i<conts.size(); ++i){
		if (IsClosed(conts[i])) throw std::runtime_error(
			"Closed contours are not allowed for "
			"contour-from-subcontours constructor");
	}
	vector<bool> fix(conts.size(), false);
	for (auto f: fixed) fix[f] = true;
	//collect endpoints
	Point f0=*First(conts[0]), f1=*First(conts[1]),
	      l0=*Last(conts[0]), l1=*Last(conts[1]);
	vector<Point> endpoints {f0, l0};
	vector<bool> reverted {false};
	double d00 = Point::meas(f0, f1);
	double d10 = Point::meas(l0, f1);
	double d01 = Point::meas(f0, l1);
	double d11 = Point::meas(l0, l1);
	if (std::min(d10, d11) > std::min(d00, d01)){
		reverted[0] = true;
		std::swap(endpoints[0], endpoints[1]);
	}
	for (int i=1; i<conts.size(); ++i){
		auto& c=conts[i];
		shared_ptr<Vertex> p1 = First(c), p2 = Last(c);
		double d0 = Point::meas(endpoints.back(), *p1); 
		double d1 = Point::meas(endpoints.back(), *p2); 
		reverted.push_back(d0 > d1);
		if (reverted.back()) std::swap(p1, p2);
		endpoints.push_back(*p2);
		auto& pcur = endpoints.back();
		auto& pprev = endpoints.end()[-2];
		if (!fix[i] && !fix[i-1]){
			if (fullshift) pcur += (pprev - *p1);
			else pprev.set(Point::Weigh(pprev, *p1, 0.5));
		} else if (fix[i] && !fix[i-1]){
			pprev.set(*p1); pcur.set(*p2);
		} else if (!fix[i] && fix[i-1]){
			if (fullshift) pcur += (pprev - *p1);
		} else if (fix[i] && fix[i-1]){
			if (pprev != *p1) throw std::runtime_error("Unconnected fixed edges detected");
		}
	}
	if (last_close){
		auto &pf = endpoints[0], &pl = endpoints.back();
		Point p = pl;
		if (!fix.back() && !fix[0]){
			p = (pf + pl) / 2;
		} else if (fix.back() && !fix[0]){
			p = pl;
		} else if (!fix.back() && fix[0]){
			p = pf;
		} else if (fix.back() && fix[0]){
			if (pf != pl) throw std::runtime_error("Unconnected fixed edges detected");
		}
		pf.set(p); pl.set(p);
	}
	//assembling
	EdgeData ret;
	for (int i=0; i<conts.size(); ++i){
		EdgeData acont;
		DeepCopy(conts[i], acont);

		auto ap = OrderedPoints(acont);
		vector<double> rw = EWeights(acont);
		int i1=i, i2=i+1;
		if (reverted[i]) std::swap(i1, i2);
		Vect v1 = endpoints[i1] - *ap[0];
		Vect v2 = endpoints[i2] - *ap.back();
		for (int j=0; j<ap.size(); ++j){
			Point p1 = *ap[j] + v1;
			Point p2 = *ap[j] + v2;
			ap[j]->set(Point::Weigh(p1, p2, rw[j]));
		}
		ret.insert(ret.end(), acont.begin(), acont.end());
	}
	ECol::Algos::MergePoints(ret);
	return Contour::Assembler::Contour1(ret, ret[0]->first().get());
}

namespace{
struct SepAssembler{
	std::list<EdgeData*> itset;
	std::map<Point*, vector<Point*>> ppmap;
	bool has_conts(){ return itset.size() > 0; }

	SepAssembler(std::list<EdgeData>& clist){
		for (auto it=clist.begin(); it!=clist.end(); ++it){
			itset.push_back(&(*it));
			auto er = ppmap.emplace(First(*it).get(), vector<Point*>());
			er.first->second.push_back((*it)[0]->last().get());
		}
		//sort ppmap by angle
		for (auto& it: ppmap){
			std::map<Point*, double> anmap;
			for (auto& it2: it.second){
				Vect v = *it2 - *it.first;
				anmap[it2] = atan2(v.y, v.x);
			}
			std::sort(it.second.begin(), it.second.end(),
					[&anmap](Point* p1, Point* p2){ return anmap[p1]<anmap[p2];} );
		}
	}
	EdgeData* take_this(EdgeData* it){
		itset.remove(it);
		return it;
	}
	EdgeData* take_any(){
		return take_this(*itset.begin());
	}
	EdgeData* take_next(EdgeData* it){
		auto fnd = ppmap.find(Last(*it).get());
		assert(fnd != ppmap.end());
		Point* pprev= it->back()->first().get();
		auto fnd2 = std::find(fnd->second.begin(), fnd->second.end(), pprev);
		assert(fnd2 != fnd->second.end());
		int ind2 = fnd2 - fnd->second.begin();
		if (ind2 == 0) ind2 = fnd->second.size()-1;
		else --ind2;
		Point* p2 = fnd->second[ind2];
		//find contour which first edge is [*it->last, p2]
		for (auto it1: itset){
			if ((*it1)[0]->first() == Last(*it) && (*it1)[0]->last().get() == p2)
				return take_this(it1);
		}
		throw std::runtime_error("failed to assemble closed contour");
	}
};
}

vector<EdgeData> cns::ExtendedSeparate(const EdgeData& _ecol){
	EdgeData ecol;
	DeepCopy(_ecol, ecol);
	//assembling subcontours
	EdgeData ecol2 = ECol::Algos::NoCrosses(ecol);
	std::vector<EdgeData> vconts = Contour::Assembler::SimpleContours(ecol2);
	std::list<EdgeData> conts;
	for (auto& c: vconts) conts.push_back(std::move(c));
	//nodes building
	std::map<Point*, int> nmap;
	for (auto it: conts) if (IsOpen(it)){
		auto er1 = nmap.emplace(First(it).get(), 0);
		auto er2 = nmap.emplace(Last(it).get(), 0);
		er1.first->second += 1;
		er2.first->second += 1;
	}
	vector<Point*> nodes;
	for (auto& it: nmap) if (it.second > 2) nodes.push_back(it.first);
	//get rid of open contours and fully closed contours
	vector<EdgeData> ret;
	for (auto it=conts.begin(); it!=conts.end();){
		if (IsClosed(*it) ||
		    std::find(nodes.begin(), nodes.end(), First(*it).get()) == nodes.end() ||
		    std::find(nodes.begin(), nodes.end(), Last(*it).get()) == nodes.end()){
			ret.push_back(EdgeData(std::move(*it)));
			it = conts.erase(it);
		} else ++it;
	}
	//double conts
	{
		int isz = conts.size();
		auto it = conts.begin();
		for (int i=0; i<isz; ++i, ++it){
			R::ReallyDirect::Permanent(*it);
			conts.push_back(EdgeData());
			DeepCopy(*it, conts.back(), 0);
			R::ReallyRevert::Permanent(conts.back());
		}
	}
	//initialize builder
	SepAssembler sep_builder(conts);
	//start assembling contours
	while (sep_builder.has_conts()){
		EdgeData bf;
		auto it = sep_builder.take_any();
		Point* p0 = First(*it).get();
		while (1){
			Algos::Connect(bf, *it);
			if (Last(*it).get() == p0) break;
			else it = sep_builder.take_next(it);
		};
		if (Contour::Area(bf) > 0){
			ret.push_back(EdgeData(std::move(bf)));
		}
	}

	return ret;
}

