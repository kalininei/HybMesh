#include "contclipping.hpp"
#include "clipper_core.hpp"

namespace ci = HMCont2D::Clip;
using namespace ci;

//two contours. direction is not taken into account
TRet ci::Intersection(const ECont& c1, const ECont& c2){
	vector<Impl::ClipperPath> p1{c1}, p2{c2};
	return Impl::ClipperPath::Intersect(p1, p2, false, false);
}

TRet ci::Union(const ECont& c1, const ECont& c2){
	vector<Impl::ClipperPath> p1{c1}, p2{c2};
	return Impl::ClipperPath::Union(p1, p2, false, false);
}

TRet ci::Difference(const ECont& c1, const ECont& c2){
	vector<Impl::ClipperPath> p1{c1}, p2{c2};
	return Impl::ClipperPath::Substruct(p1, p2, false, false);
}


//tree and contour: contours direction is not taken into account
TRet ci::Intersection(const ETree& c1, const ECont& c2){
	_THROW_NOT_IMP_;
}

TRet ci::Union(const ETree& c1, const ECont& c2){
	vector<Impl::ClipperPath> p1, p2{c2};
	for (auto& n: c1.nodes){
		p1.push_back(Impl::ClipperPath(*n));
	}
	return Impl::ClipperPath::Union(p1, p2, true, false);
}

TRet ci::Difference(const ETree& c1, const ECont& c2){
	vector<Impl::ClipperPath> p1, p2{c2};
	for (auto& n: c1.nodes){
		p1.push_back(Impl::ClipperPath(*n));
	}
	return Impl::ClipperPath::Substruct(p1, p2, true, false);
}

TRet ci::Difference(const ECont& c1, const ETree& c2){
	_THROW_NOT_IMP_;
}


//two trees
TRet ci::Intersection(const ETree& c1, const ETree& c2){
	_THROW_NOT_IMP_;
}

TRet ci::Union(const ETree& c1, const ETree& c2){
	_THROW_NOT_IMP_;
}

TRet ci::Difference(const ETree& c1, const ETree& c2){
	_THROW_NOT_IMP_;
}


//multiple contours operation.
//Direction of each contour is not taken into account.
TRet ci::Intersection(const vector<ECont>& cont){
	_THROW_NOT_IMP_;
}

TRet ci::Union(const vector<ECont>& cont){
	if (cont.size() == 0) return TRet();
	vector<Impl::ClipperPath> p1, zero;
	for (auto& c: cont) p1.push_back(Impl::ClipperPath(c));
	return Impl::ClipperPath::Union(p1, zero, false, false);
}

TRet ci::Difference(const ETree& c1, const vector<ECont>& cont){
	if (cont.size() == 0) return TRet::DeepCopy(c1);
	vector<Impl::ClipperPath> p1, p2;
	for (auto& c: c1.nodes) p1.push_back(Impl::ClipperPath(*c));
	for (auto& c: cont) p2.push_back(Impl::ClipperPath(c));
	return Impl::ClipperPath::Substruct(p1, p2, true, false);
}

void ci::Heal(TRet& c1){
	//make union to connect adjacent areas if they have appeared
	if (c1.cont_count() > 1){
		vector<Impl::ClipperPath> p1, zero;
		for (auto& c: c1.nodes) p1.push_back(Impl::ClipperPath(*c));
		c1 = Impl::ClipperPath::Union(p1, zero, true, false);
	}
	//remove zero length edges
	for (auto& path: c1.nodes){
		auto pts = path->corner_points();
		if (pts.size() < 3){ 
			//invalid contour. Will be removed.
			path->data.resize(2);
			continue;
		}
		//doubled points
		std::set<int> not_needed;
		for (int i=0; i<pts.size(); ++i){
			auto iprev = (i==0) ? pts.size()-1 : i-1;
			if (Point::meas(*pts[iprev], *pts[i])<sqr(geps)) not_needed.insert(i);
		}
		aa::remove_entries(pts, not_needed);
		not_needed.clear();
		//collinear nodes
		for (int i=0; i<pts.size(); ++i){
			auto iprev = (i==0) ? pts.size()-1 : i-1;
			auto inext = (i==pts.size()-1) ? 0 : i+1;
			if (fabs(triarea(*pts[iprev], *pts[i], *pts[inext]))<geps){
				not_needed.insert(i);
			}
		}
		aa::remove_entries(pts, not_needed);
		not_needed.clear();
		//wright simplified contour
		auto c = HMCont2D::Constructor::ContourFromPoints(pts, true);
		path->data = c.data;
	}
	//remove invalid contours
	std::set<int> not_needed;
	for (int i=0; i<c1.cont_count(); ++i)
		if (c1.nodes[i]->size()<3) not_needed.insert(i);
	aa::remove_entries(c1.nodes, not_needed);
	//updates topology
	c1.UpdateTopology();
	//refill data
	c1.data.clear();
	for (auto& n: c1.nodes) std::copy(n->begin(), n->end(), std::back_inserter(c1.data));
}


