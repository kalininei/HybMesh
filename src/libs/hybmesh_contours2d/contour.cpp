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
	auto pnt = ordered_points();
	vector<Point*> ret;
	for (int i=0; i<pnt.size()-1; ++i){
		Point *pprev, *p, *pnext;
		if (i == 0){
			if (is_closed()) pprev = pnt[pnt.size() - 2];
			else continue;
		} else pprev = pnt[i-1];
		p = pnt[i];
		pnext = pnt[i+1];
		//build triangle. if its area = 0 -> that is not a corner point
		if (fabs(triarea(*pprev, *p, *pnext)) > geps)
			ret.push_back(p);
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

std::tuple<bool, Point*>
Contour::GuaranteePoint(const Point& p, PCollection& pcol){
	_DUMMY_FUN_;
	std::tuple<bool, Point*> ret;
	std::get<0>(ret) = false;
	std::get<1>(ret) = FindClosestNode(*this, p);
	return ret;
}


Container<ContourTree> Contour::Offset(const Contour& source, double delta, OffsetTp tp){
	Impl::ClipperPath cp(source);
	return cp.Offset(delta, tp);
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

Contour Contour::Assemble(const Contour& con, Point* pnt_start, Point* pnt_end){
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

Contour Contour::Assemble(const ECollection& col, Point* pnt_start, Point* pnt_end){
	//1) get contour which contains pnt_start
	Contour x = Contour::Assemble(col, pnt_start);
	//2) cut it
	return Contour::Assemble(x, pnt_start, pnt_end);
}

Contour Contour::Assemble(const ECollection& col, Point* pnt_start){
	std::list<shared_ptr<Edge>> eset(col.begin(), col.end());
	
	//finds and removes from eset Edge with point
	auto fnded = [&eset](Point* p) -> shared_ptr<Edge>{ 
		auto fnd = std::find_if(eset.begin(), eset.end(),
			[&p](shared_ptr<Edge> e){ return e->contains(p); });
		if (fnd == eset.end()) return nullptr;
		auto ret = *fnd;
		eset.erase(fnd);
		return ret;
	};

	//finds another point of edge
	auto nextnd = [](Edge* e, Point* p)->Point*{ return (p == e->pstart) ? e->pend : e->pstart; };

	//builds edges in any direction
	auto build_half_cont = [&](Point* p1){
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


PCollection Contour::WeightPoints(const Contour& p, vector<double> vw){
	PCollection ret;
	auto pseq = p.ordered_points();
	std::sort(vw.begin(), vw.end());
	assert(!ISLOWER(vw[0], 0) && !ISGREATER(vw.back(), 1.0));

	vector<double> elens {0};
	for (auto& e: p.data) elens.push_back(Point::dist(*e->pstart, *e->pend));
	std::partial_sum(elens.begin(), elens.end(), elens.begin());
	double len = elens.back();
	std::transform(vw.begin(), vw.end(), vw.begin(), [&len](double x){ return x*len; });

	int icur=1;
	for (auto w: vw){
		while (elens[icur]<w && icur<pseq.size()) ++icur;
		Point *pprev = pseq[icur-1], *pnext = pseq[icur];
		double t = (w - elens[icur-1])/(elens[icur] - elens[icur-1]);
		ret.add_value( (*pseq[icur-1]) * (1-t) + (*pseq[icur]) * t );
	}
	return ret;
}

