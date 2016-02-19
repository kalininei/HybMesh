#include "contour.hpp"
#include "tree.hpp"
#include "clipper_core.hpp"
#include "debug_cont2d.hpp"

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

vector<Point*> Contour::corner_points() const{
	bool cl = is_closed();
	auto pnt = ordered_points();
	vector<Point*> ret;
	if (!cl) ret.push_back(pnt[0]);
	for (int i=0; i<pnt.size()-1; ++i){
		Point *pprev, *p, *pnext;
		if (i == 0){
			if (cl) pprev = pnt[pnt.size() - 2];
			else continue;
		} else pprev = pnt[i-1];
		p = pnt[i];
		pnext = pnt[i+1];
		//build triangle. if its area = 0 -> that is not a corner point
		//if (fabs(triarea(*pprev, *p, *pnext)) > geps)
			//ret.push_back(p);
		double ksi = Point::meas_section(*p, *pprev, *pnext);
		if (ksi>=geps*geps){
			ret.push_back(p);
		}
	}
	if (!cl) ret.push_back(pnt.back());
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
	//e1
	HMCont2D::Edge *e1 = data[i].get(), *e2;
	//e2
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
		if (e1->pstart == e2->pstart) return {e1->pend, e1->pstart, e2->pend};
		if (e1->pstart == e2->pend) return {e1->pend, e1->pstart, e2->pstart};
		if (e1->pend == e2->pstart) return {e1->pstart, e1->pend, e2->pend};
		return {e1->pstart, e1->pend, e2->pstart};
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
				a.pnext = pts[2];
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

Point Contour::ClosestPoint(const Point& p) const{
	auto ce = FindClosestEdge(*this, p);
	if (ISEQ(std::get<2>(ce), 0)){
		return *std::get<0>(ce)->pstart;
	} else if (ISEQ(std::get<2>(ce), 1)){
		return *std::get<0>(ce)->pend;
	} else {
		Point* p1 = std::get<0>(ce)->pstart;
		Point* p2 = std::get<0>(ce)->pend;
		return Point::Weigh(*p1, *p2, std::get<2>(ce));
	}
}

Container<ContourTree> Contour::Offset(const Contour& source, double delta, OffsetTp tp){
	Impl::ClipperPath cp(source);
	if (source.is_closed() && Contour::Area(source) < 0) delta = -delta;
	return cp.Offset(delta, tp);
};

Container<Contour> Contour::Offset1(const Contour& source, double delta){
	Container<ContourTree> ans;
	if (source.is_closed()) ans = Offset(source, delta, OffsetTp::CLOSED_POLY);
	else ans = Offset(source, delta, OffsetTp::OPEN_ROUND);
	assert(ans.cont_count() == 1);
	return Container<Contour>::DeepCopy(*ans.get_contour(0));
};

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

namespace{
Vect smoothed_direction_core(const vector<Point*>& p){
	assert(p.size() > 1);
	Vect ret1 = *p[1] - *p[0];
	Vect ret2 = *p.back() - *p[0];
	double a = Angle(ret1, Point(0, 0), ret2);
	Vect ret;
	if (a<M_PI/4 || a>3*M_PI/4) ret = ret1;
	else ret = ret2;
	vecNormalize(ret);
	return ret;
}

}

Vect Contour::SmoothedDirection(const Contour& c, Point* p, int direction, double len){
	std::list<Point> apoints; //storage for additional points
	vector<Point*> chosen = c.ordered_points();
	//1. place p into chosen array
	int pindex = -1;
	//try to find exact match
	auto fnd = std::find(chosen.begin(), chosen.end(), p);
	if (fnd == chosen.end()) fnd = std::find_if(chosen.begin(), chosen.end(),
			[&p](Point* p1){ return *p1 == *p; });
	if (fnd != chosen.end()) pindex = fnd - chosen.begin();
	//if failed -> place p on nearest edge
	if (pindex < 0){
		auto ce = FindClosestEdge(c, *p);
		Edge* e = std::get<0>(ce);
		double &w = std::get<2>(ce);
		if (ISEQ(w, 0)) return SmoothedDirection(c,  e->pstart, direction, len);
		if (ISEQ(w, 1)) return SmoothedDirection(c,  e->pend, direction, len);
		apoints.push_back(Point::Weigh(*e->pstart, *e->pend, w));
		pindex = c.get_index(e) + 1;
		chosen.insert(chosen.begin() + pindex, &apoints.back());
	}
	
	//2. revert according to direction
	if (direction == -1){
		std::reverse(chosen.begin(), chosen.end());
		pindex = chosen.size()- 1 - pindex;
	}
	//3. for closed contours set p as the first point
	//and guarantee len > 0.25 of full length. Hence we can remove edge before pindex.
	if (c.is_closed()){
		chosen.pop_back();
		for (int i=0; i<pindex; ++i) chosen.push_back(chosen[i]);
		chosen = vector<Point*>(chosen.begin()+pindex, chosen.end());
		pindex = 0;
		len = std::min(len, 0.25*c.length());
	} else {
	//4. if full length of open contour is less then len treat all contour
		if (len>=c.length()) return smoothed_direction_core(chosen);
	}
	//5. go forward till we can
	int iend = pindex;
	double usedx = 0, last_len=0;
	while (iend<chosen.size()-1 && ISGREATER(len, usedx)){
		++iend;
		last_len = Point::dist(*chosen[iend-1], *chosen[iend]);
		usedx+=last_len;
	}
	//if we've gone too far change chosen[iend] to weighted point to provide len.
	if (ISGREATER(usedx, len)){
		double w = (usedx - len)/last_len;
		apoints.push_back(Point::Weigh(*chosen[iend-1], *chosen[iend], w));
		chosen[iend] = &apoints.back();
	}
	//6 .go backward if necessary
	int istart = pindex;
	if (ISGREATER(len, usedx)){
		//_THROW_NOT_IMP_;
	}
	
	//7. leave only [istart-iend] points and calculate
	chosen = vector<Point*>(chosen.begin()+istart,  chosen.begin() + iend+1);
	return smoothed_direction_core(chosen);
}


Contour Contour::Partition(double step, const Contour& contour, PCollection& pstore, PartitionTp tp){
	switch (tp){
		case PartitionTp::IGNORE_ALL:
			return Partition(step, contour, pstore);
		case PartitionTp::KEEP_ALL:
			return Partition(step, contour, pstore, contour.all_points());
		case PartitionTp::KEEP_SHAPE:
			return Partition(step, contour, pstore, contour.corner_points());
	};
}

namespace {

//build a new contour based on input contour begin/end points
Contour partition_core(double step, const Contour& contour, PCollection& pstore){
	double len = contour.length();
	if (len<1.5*step){
		if (!contour.is_closed()){
			Contour ret;
			ret.add_value(Edge{contour.first(), contour.last()});
			return ret;
		} else {
			step = len/3.1;
		}
	}
	//calculates new points
	int n = (int)round(len/step);
	vector<double> w;
	for (int i=1; i<n; ++i) w.push_back((double)i/n);
	PCollection wp = Contour::WeightPoints(contour, w);
	//add new poins to pstore
	pstore.Unite(wp);
	//Construct new contour
	vector<Point*> cpoints;
	cpoints.push_back(contour.first());
	std::transform(wp.begin(), wp.end(), std::back_inserter(cpoints),
			[](shared_ptr<Point> x){ return x.get(); });
	cpoints.push_back(contour.last());
	return Constructor::ContourFromPoints(cpoints);
}

//returns repartitioned copy of a contour 
Contour partition_core(double step, const Contour& contour, PCollection& pstore, const std::list<Point*>& keep){
	auto it0 = keep.begin(), it1 = std::next(it0);
	Contour ret;
	while (it1 != keep.end()){
		Contour sub = Contour::Assemble(contour, *it0, *it1);
		Contour psub = partition_core(step, sub, pstore);
		ret.Unite(psub);
		//if sub.size() == 1 then its direction is not defined
		//so we need to check the resulting direction.
		//We did it after second Unition when direction matters
		if (it0 == std::next(keep.begin())){
			if (ret.last() != *it1) ret.Reverse();
		}

		++it0; ++it1;
	}
	return ret;
}

}

Contour Contour::Partition(double step, const Contour& contour, PCollection& pstore,
		const std::vector<Point*>& keepit){
	//sort points in keepit and add start, end point there
	vector<Point*> orig_pnt = contour.ordered_points();
	std::set<Point*> keepset(keepit.begin(), keepit.end());

	std::list<Point*> keep_sorted;
	std::copy_if(orig_pnt.begin(), orig_pnt.end(), std::back_inserter(keep_sorted),
		[&keepset](Point* p){ return keepset.find(p) != keepset.end(); }
	);
	
	//add first, last points
	if (!contour.is_closed()){
		if (*keep_sorted.begin() != orig_pnt[0])
			keep_sorted.push_front(orig_pnt[0]);
		if (*keep_sorted.rbegin() != orig_pnt.back())
			keep_sorted.push_back(orig_pnt.back());
	} else if (keep_sorted.size() == 1){
		keep_sorted.push_back(*keep_sorted.begin());
	} else if (keep_sorted.size() == 0){
		keep_sorted.push_back(orig_pnt[0]);
		keep_sorted.push_back(orig_pnt[0]);
	}
	
	//call core procedure
	return partition_core(step, contour, pstore, keep_sorted);
}

bool Contour::DoIntersect(const Contour& c1, const Contour& c2){
	auto bbox1 = Contour::BBox(c1);
	auto bbox2 = Contour::BBox(c2);
	if (!bbox1.has_common_points(bbox2)) return false;
	auto c = HMCont2D::Clip::Union(c1, c2);
	if (c.cont_count() > 1) return false;
	else return true;
}

bool Contour::DoReallyIntersect(const Contour& c1, const Contour& c2){
	auto bbox1 = Contour::BBox(c1);
	auto bbox2 = Contour::BBox(c2);
	if (!bbox1.has_common_points(bbox2)) return false;
	auto c = HMCont2D::Clip::Union(c1, c2);
	if (c.cont_count() > 1) return false;
	double a1 = fabs(HMCont2D::Area(c1));
	double a2 = fabs(HMCont2D::Area(c2));
	double suma = fabs(HMCont2D::Area(*c.nodes[0]));
	return (!ISZERO(a1+a2 - suma));
}

namespace{

Contour assemble_core(const Contour& con, int estart, int eend){
	Contour ret;
	for (int i=estart; i<eend; ++i){
		int k = i<con.size() ? i : i-con.size();
		ret.add_value(con.data[k]);
	}
	return ret;
}

}

Contour Contour::Assemble(const Contour& con, const Point* pnt_start, const Point* pnt_end){
	//find indicies of pnt_start/pnt_end
	auto op = con.ordered_points();
	int i0 = std::find_if(op.begin(), op.end(), [&pnt_start](Point* p){ return p == pnt_start; })
		- op.begin();
	int i1 = std::find_if(op.begin(), op.end(), [&pnt_end](Point* p){ return p == pnt_end; })
		- op.begin();
	assert(i0 < op.size() && i1 < op.size());

	if (con.is_closed()){
		if (i0 == i1) i1 = i0 + con.size();
		if (i1 == 0) i1 = con.size();
	}
	
	if (i1 >= i0) return assemble_core(con, i0, i1);
	else {
		if (con.is_closed()) return assemble_core(con, i0, i1 + con.size());
		else{
			//if reversed non-closed
			Contour r = assemble_core(con, i1, i0);
			r.Reverse();
			return r;
		}
	}
}

Contour Contour::Assemble(const ECollection& col, const Point* pnt_start, const Point* pnt_end){
	//1) get contour which contains pnt_start
	Contour x = Contour::Assemble(col, pnt_start);
	//2) cut it
	return Contour::Assemble(x, pnt_start, pnt_end);
}

Contour Contour::Assemble(const ECollection& col, const Point* pnt_start){
	std::list<shared_ptr<Edge>> eset(col.begin(), col.end());
	
	//finds and removes from eset Edge with point
	auto fnded = [&eset](const Point* p) -> shared_ptr<Edge>{ 
		auto fnd = std::find_if(eset.begin(), eset.end(),
			[&p](shared_ptr<Edge> e){ return e->contains(p); });
		if (fnd == eset.end()) return nullptr;
		auto ret = *fnd;
		eset.erase(fnd);
		return ret;
	};

	//finds another point of edge
	auto nextnd = [](Edge* e, const Point* p)->Point*{ return (p == e->pstart) ? e->pend : e->pstart; };

	//builds edges in any direction
	auto build_half_cont = [&](const Point* p1){
		Contour r;
		while(1){
			auto a = fnded(p1);
			if (a) {
				r.add_value(a);
				p1 = nextnd(a.get(), p1);
			} else break;
		}
		return r;
	};

	auto c1 = build_half_cont(pnt_start);
	auto c2 = build_half_cont(pnt_start);
	c1.Unite(c2);
	return c1;
}

Contour Contour::Assemble(const Contour& col, const Point* pnt_start, int direction, double len){
	assert(col.contains_point(pnt_start));
	if (col.size() == 1) return Contour(col);
	if (direction == -1){
		Contour col2(col);
		col2.Reverse();
		return Assemble(col2, pnt_start, 1, len);
	}

	//assemble ordered points taking into account the direction
	vector<Point*> op = col.ordered_points();
	int index = std::find(op.begin(), op.end(), pnt_start) - op.begin();
	//if is_closed place start point at the begining
	if (col.is_closed()){
		op.pop_back();
		std::rotate(op.begin(), op.begin() + index, op.end());
		index = 0;
		op.push_back(op[0]);
	}
	//go from start till distance is more then len in order to find
	//last point
	double dist = 0;
	Point* pnt_end;
	do{
		pnt_end = op[++index];
		dist += Point::dist(*op[index-1], *op[index]);
		if (ISEQGREATER(dist, len)) break;
	} while (index<op.size());

	return Assemble(col, pnt_start, pnt_end);
}

Contour Contour::Simplified(const Contour& cont){
	auto p = cont.corner_points();
	return HMCont2D::Constructor::ContourFromPoints(p, cont.is_closed());
}

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

namespace{
vector<std::tuple<bool, Point, double, double>>
cross_core(const Contour& c1, const Contour& c2, bool is1){
	vector<std::tuple<bool, Point, double, double>> retv;
	vector<Point*> op1 = c1.ordered_points();
	vector<Point*> op2 = c2.ordered_points();

	auto lens1 = ECollection::ELengths(c1), lens2 = ECollection::ELengths(c2);
	double flen1 = std::accumulate(lens1.begin(), lens1.end(), 0.0);
	double flen2 = std::accumulate(lens2.begin(), lens2.end(), 0.0);

	auto addcross = [&](Point p, double w1, double w2){
		retv.push_back(std::make_tuple(true, p, w1, w2));
	};
	double ksieta[2];
	double L1=0.0, L2=0.0;
	for (int i=0; i<op1.size()-1; ++i){
		L2 = 0;
		for (int j=0; j<op2.size()-1; ++j){
			if (SectCross(*op1[i], *op1[i+1], *op2[j], *op2[j+1], ksieta)){
				addcross(
					Point::Weigh(*op1[i], *op1[i+1], ksieta[0]),
					(L1 + lens1[i]*ksieta[0])/flen1,
					(L2 + lens2[j]*ksieta[1])/flen2
				);
				if (is1) return retv;
			}
			L2+=lens2[j];
		}
		L1 += lens1[i];
	}
	if (retv.size() < 2) return retv;

	//clear doublicates
	TCoordSet<double> w1;
	vector<std::tuple<bool, Point, double, double>> ret;
	for (auto& v: retv){
		if (w1.find(std::get<2>(v)) == w1.end()){
			w1.insert(std::get<2>(v));
			ret.push_back(v);
		}
	}
	if (ret.size() == 1) return ret;

	return ret;
}
}

std::tuple<bool, Point, double, double>
Contour::Cross(const Contour& c1, const Contour& c2){
	auto retv = cross_core(c1, c2, true);
	if (retv.size() == 0) return std::make_tuple(false, Point(0,0), 0.0, 0.0);
	else return retv[0];
}

vector<std::tuple<bool, Point, double, double>>
Contour::CrossAll(const Contour& c1, const Contour& c2){
	return cross_core(c1, c2, false);
}

Container<Contour> Contour::CutByWeight(const Contour& source, double w1, double w2){
	Container<Contour> tmp;
	Container<Contour>::DeepCopy(source, tmp);
	PCollection wp = WeightPoints(tmp, {w1, w2});
	auto p1 = tmp.GuaranteePoint(*wp.point(0), tmp.pdata);
	auto p2 = tmp.GuaranteePoint(*wp.point(1), tmp.pdata);
	Contour rc = Contour::Assemble(tmp, std::get<1>(p1), std::get<1>(p2));
	Container<Contour> ret;
	Container<Contour>::DeepCopy(rc, ret);
	return ret;
}

Container<Contour> Contour::CutByLen(const Contour& source, double len1, double len2){
	double len = source.length();
	return CutByWeight(source, len1/len, len2/len);
}






