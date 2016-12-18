#include "constructor.hpp"
#include "edges.hpp"
#include "spmat.hpp"
#include "contour.hpp"

using namespace HMCont2D;
namespace cns = Constructor;

Container<Contour> cns::Circle(int N, double rad, Point cnt){
	Point poc = cnt + Point(rad, 0);
	return Circle(N, cnt, poc);
};

HMCont2D::Container<Contour> cns::Circle(int N, Point cnt, Point poc){
	typedef Container<Contour> TRes;
	TRes res;
	Point p = poc-cnt;
	double rad = Point::dist(poc, cnt);
	double angle0 = atan2(p.y, p.x);
	for (int i=0; i<N; ++i){
		double angle = 2*M_PI/N*i + angle0;
		Point a(rad*cos(angle), rad*sin(angle));
		a += cnt;

		res.pdata.add_value(a);

		if (i>0) res.add_value(Edge(res.point(i-1), res.point(i)));
	}
	res.add_value(Edge(res.point(N-1), res.point(0)));
	return res;
}

Contour cns::ContourFromPoints(const vector<Point*>& pnt, bool force_closed){
	Contour ret;
	if (pnt.size()<2) return ret;
	for (int i=0; i<(int)pnt.size() - 1; ++i){
		ret.add_value(Edge(pnt[i], pnt[i+1]));
	}
	if (force_closed && pnt[0]!=pnt.back()){
		ret.add_value(Edge(pnt.back(), pnt[0]));
	}
	return ret;
}

Container<Contour> cns::ContourFromPoints(vector<double> pnt, bool force_closed){
	vector<Point> p;
	for (int i=0; i<(int)pnt.size(); i+=2) p.push_back(Point(pnt[i], pnt[i+1]));
	return cns::ContourFromPoints(p, force_closed);
}

Container<Contour> cns::ContourFromPoints(vector<Point> pnt, bool force_closed){
	PCollection pc;
	for (auto& it: pnt) pc.add_value(it);
	vector<Point*> p = pc.pvalues();
	if (pnt[0] == pnt.back()){
		pc.pop_value();
		p.back() = p[0];
	}
	Contour cc = cns::ContourFromPoints(p, force_closed);
	Container<Contour> ret;
	ret.data = cc.data;
	ret.pdata = pc;
	return ret;	
}

Contour cns::ContourFromPoints(const HMCont2D::PCollection& dt, bool force_closed){
	return cns::ContourFromPoints(dt.pvalues(), force_closed);
}

Container<Contour> cns::ContourFromBBox(BoundingBox bbox){
	return ContourFromPoints(bbox.FourPoints(), true);
}

Container<Contour> cns::CutContour(const HMCont2D::Contour& cont,
		const Point& pstart, int direction, double len){
	Container<Contour> ret;
	//if pstart lies on cont -> use simple assembling
	vector<Point*> op = cont.ordered_points();
	auto fnd = std::find_if(op.begin(), op.end(), [&pstart](Point* p){return pstart == *p;});
	if (fnd != op.end()){
		Contour c2 = Assembler::Contour1(cont, *fnd, direction, len);
		Container<Contour>::DeepCopy(c2, ret);
		if (*ret.first() != pstart) ret.ReallyReverse();
		return ret;
	}
	//if not-> place point and try once again
	Container<Contour> c2 = Container<Contour>::DeepCopy(cont);
	c2.GuaranteePoint(pstart, c2.pdata);
	return cns::CutContour(c2, pstart, direction, len);
}

Container<Contour> cns::CutContour(const HMCont2D::Contour& cont,
		const Point& pstart, const Point& pend){
	PCollection tmpp;
	Contour tmpc = cont;
	auto p1 = std::get<1>(tmpc.GuaranteePoint(pstart, tmpp));
	auto p2 = std::get<1>(tmpc.GuaranteePoint(pend, tmpp));
	tmpc = Assembler::Contour1(tmpc, p1, p2);
	Container<Contour> ret;
	Container<Contour>::DeepCopy(tmpc, ret);
	if (*ret.first() != *p1) ret.ReallyReverse();
	return ret;
}

Container<Contour> cns::CutContourByWeight(const Contour& source, double w1, double w2){
	PCollection wp = Contour::WeightPoints(source, {w1, w2});
	return CutContour(source, wp.value(0), wp.value(1));
}
Container<Contour> cns::CutContourByLen(const Contour& source, double len1, double len2){
	double len = source.length();
	return CutContourByWeight(source, len1/len, len2/len);
}

HMCont2D::Container<HMCont2D::ECollection> cns::ECol(const vector<Point>& pnt, const vector<int>& eds){
	HMCont2D::Container<HMCont2D::ECollection> ret;
	ret.pdata.data.reserve(pnt.size());
	ret.data.reserve(eds.size()/2);
	for (auto& p: pnt) ret.pdata.add_value(p);

	for (int i=0; i<eds.size()/2; ++i){
		int i1 = eds[2*i];
		int i2 = eds[2*i+1];
		ret.add_value(HMCont2D::Edge(ret.point(i1), ret.point(i2)));
	}
	return ret;
}
HMCont2D::Container<HMCont2D::ECollection> cns::ECol(int npnt, int neds, double* pnt, int* eds){
	vector<int> evec(eds, eds + 2*neds);
	vector<Point> pvec; pvec.reserve(npnt);
	for (int i=0; i<npnt; ++i){
		pvec.push_back(Point(pnt[2*i], pnt[2*i+1]));
	}
	return ECol(pvec, evec);
}

HMCont2D::Container<HMCont2D::Contour> cns::PerturbedContour(Point p1, Point p2, int npart,
		std::function<double(double)> perturbation){
	vector<Point> pp(npart+1);
	pp[0] = p1;
	pp.back() = p2;
	Vect pvec = p2 - p1;
	for (int i=1; i<npart; ++i){
		double t = (double)i/(npart);
		double m = perturbation(t);
		Vect perpvec(-pvec.y, pvec.x);
		if (m<0) perpvec *= -1.0;
		vecSetLen(perpvec, fabs(m));
		pp[i] = Point::Weigh(p1, p2, t);
		pp[i] += perpvec;
	}
	return cns::ContourFromPoints(pp);
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
HMCont2D::Contour cns::Spline(const vector<Point*>& pnt, HMCont2D::PCollection& pcol, int nedges, bool force_closed){
	//force closed option
	if (force_closed && pnt[0] != pnt.back()){
		vector<Point*> pnew = pnt;
		if (*pnt[0] != *pnt.back()) pnew.push_back(pnt[0]);
		else pnew.back() = pnt[0];
		return Spline(pnew, pcol, nedges, false);
	}

	int nseg = pnt.size() - 1;
	vector<double> double_nsum;
	vector<double> vlen(1, 0.0);
	for (int i=0; i<nseg; ++i){
		double len = Point::dist(*pnt[i+1], *pnt[i]);
		vlen.push_back(vlen.back() + len);
		double_nsum.push_back(len);
	}
	for (int i=0; i<nseg; ++i) double_nsum[i] *= (nedges/vlen.back());
	vector<int> nsum = Algos::RoundVector(double_nsum, vector<int>(nseg, 1));

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
	vector<Point*> pp;
	for (int i=0; i<nseg; ++i){
		pp.push_back(pnt[i]);
		auto& fx = xcubic[i];
		auto& fy = ycubic[i];
		for (int j=1; j<nsum[i]; ++j){
			double t = (double)j/(nsum[i]);
			double x = fx[0]*t*t*t + fx[1]*t*t + fx[2]*t + fx[3];
			double y = fy[0]*t*t*t + fy[1]*t*t + fy[2]*t + fy[3];
			shared_ptr<Point> p(new Point(x, y));
			pcol.add_value(p);
			pp.push_back(p.get());
		}
	}
	pp.push_back(pnt.back());

	return cns::ContourFromPoints(pp, false);
}

HMCont2D::Container<HMCont2D::Contour> cns::Spline(const vector<Point>& pnt, int nedges, bool force_closed){
	vector<Point*> p2;
	for (auto& p: pnt) p2.push_back(const_cast<Point*>(&p));
	if (*p2[0] == *p2.back()) p2.back() = p2[0];
	HMCont2D::PCollection pcol;
	auto cont = Spline(p2, pcol, nedges, force_closed);
	HMCont2D::Container<HMCont2D::Contour> ret;
	HMCont2D::Container<HMCont2D::Contour>::DeepCopy(cont, ret);
	return ret;
}

HMCont2D::Container<HMCont2D::Contour> cns::Spline(const vector<double>& pnt, int nedges, bool force_closed){
	vector<Point> pv;
	for (int i=0; i<pnt.size(); i+=2) pv.push_back(Point(pnt[i], pnt[i+1]));
	return Spline(pv, nedges, force_closed);
}

HMCont2D::Container<HMCont2D::Contour> cns::ContourFromContours(const vector<Contour>& conts,
		bool last_close, bool fullshift, std::set<int> fixed){
	//checks
	for (int i=0; i<conts.size(); ++i){
		if (conts[i].is_closed()) throw std::runtime_error(
			"Closed contours are not allowed for "
			"contour-from-subcontours constructor");
	}
	vector<bool> fix(conts.size(), false);
	for (auto f: fixed) fix[f] = true;
	//collect endpoints
	vector<Point> endpoints {*conts[0].first(), *conts[0].last()};
	vector<bool> reverted {false};
	double d00 = Point::meas(*conts[0].first(), *conts[1].first());
	double d10 = Point::meas(*conts[0].last(), *conts[1].first());
	double d01 = Point::meas(*conts[0].first(), *conts[1].last());
	double d11 = Point::meas(*conts[0].last(), *conts[1].last());
	if (std::min(d10, d11) > std::min(d00, d01)){
		reverted[0] = true;
		std::swap(endpoints[0], endpoints[1]);
	}
	for (int i=1; i<conts.size(); ++i){
		auto& c=conts[i];
		Point *p1 = c.first(), *p2 = c.last();
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
	HMCont2D::Container<HMCont2D::ECollection> ret;
	for (int i=0; i<conts.size(); ++i){
		HMCont2D::Container<HMCont2D::Contour> acont;
		HMCont2D::Container<HMCont2D::Contour>::DeepCopy(conts[i], acont);

		vector<Point*> ap = acont.ordered_points();
		vector<double> rw = HMCont2D::Contour::EWeights(acont);
		int i1=i, i2=i+1;
		if (reverted[i]) std::swap(i1, i2);
		Vect v1 = endpoints[i1] - *ap[0];
		Vect v2 = endpoints[i2] - *ap.back();
		for (int j=0; j<ap.size(); ++j){
			Point p1 = *ap[j] + v1;
			Point p2 = *ap[j] + v2;
			ap[j]->set(Point::Weigh(p1, p2, rw[j]));
		}
		ret.Unite(acont);
	}
	HMCont2D::Algos::MergePoints(ret);
	HMCont2D::Algos::DeleteUnusedPoints(ret);
	HMCont2D::Container<HMCont2D::Contour> ret2;
	HMCont2D::Contour rcont = HMCont2D::Assembler::Contour1(ret, ret.data[0]->pstart);
	ret2.pdata.data = std::move(ret.pdata.data);
	ret2.data = std::move(rcont.data);
	assert(ret2.check_connectivity());
	return ret2;
}
