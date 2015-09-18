#include "constructor.hpp"
#include "edges.hpp"
#include "contour.hpp"

using namespace HMCont2D;
namespace cns = Constructor;

Container<Contour> cns::Circle(int N, double rad, Point cnt){
	typedef Container<Contour> TRes;
	TRes res;
	for (int i=0; i<N; ++i){
		double angle = 2*M_PI/N*i;
		Point a(rad*cos(angle), rad*sin(angle));
		a += cnt;

		res.pdata.add_value(a);

		if (i>0) res.add_value(Edge(res.point(i-1), res.point(i)));
	}
	res.add_value(Edge(res.point(N-1), res.point(0)));
	return res;
};

Contour cns::ContourFromPoints(const vector<Point*>& pnt, bool force_closed){
	Contour ret;
	for (int i=0; i<pnt.size() - 1; ++i){
		ret.add_value(Edge(pnt[i], pnt[i+1]));
	}
	if (force_closed && pnt[0]!=pnt.back()){
		ret.add_value(Edge(pnt.back(), pnt[0]));
	}
	return ret;
}

Container<Contour> cns::ContourFromPoints(vector<double> pnt, bool force_closed){
	vector<Point> p;
	for (int i=0; i<pnt.size(); i+=2) p.push_back(Point(pnt[i], pnt[i+1]));
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

Container<Contour> cns::CutContour(const HMCont2D::Contour& cont, const Point& pstart, int direction, double len){
	Container<Contour> ret;
	//if pstart lies on cont -> use simple assembling
	vector<Point*> op = cont.ordered_points();
	auto fnd = std::find_if(op.begin(), op.end(), [&pstart](Point* p){return pstart == *p;});
	if (fnd != op.end()){
		Contour c2 = Contour::Assemble(cont, *fnd, direction, len);
		Container<Contour>::DeepCopy(c2, ret);
		return ret;
	}
	//place pstart on contour then use assembling
	Contour c2 = cont;
	auto eres = ECollection::FindClosestEdge(c2, pstart);
	Edge* ed = std::get<0>(eres);
	int index = c2.get_index(ed);
	c2.RemoveAt({index});
	Point pp = pstart;
	c2.AddAt(index, { std::make_shared<Edge>(ed->pstart, &pp), std::make_shared<Edge>(ed->pstart, &pp) });
	Container<Contour>::DeepCopy(c2, ret);
	return ret;
}

