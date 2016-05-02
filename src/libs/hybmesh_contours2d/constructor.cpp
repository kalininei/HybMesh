#include "constructor.hpp"
#include "edges.hpp"
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
