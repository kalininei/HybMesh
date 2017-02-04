#include "contclipping.hpp"
#include "clipper_core.hpp"
#include "gpc_core.hpp"
#include "cont_assembler.hpp"

namespace ci = HM2D::Contour::Clip;
using namespace ci;

namespace {
TRet assign_btypes(const ECont& c1, const ECont& c2, TRet&& ret){
	ECont ac = c1;
	ac.insert(ac.end(), c2.begin(), c2.end());
	ECont ae = ret.alledges();
	HM2D::ECol::Algos::AssignBTypes(ac, ae);
	return ret;
}
TRet assign_btypes(const vector<ECont>& c1, TRet&& ret){
	if (c1.size() == 0) return ret;
	ECont ac = c1[0];
	for (int i=1; i<c1.size(); ++i){
		ac.insert(ac.end(), c1[i].begin(), c1[i].end());
	}
	ECont ae = ret.alledges();
	HM2D::ECol::Algos::AssignBTypes(ac, ae);
	return ret;
}
TRet assign_btypes(const ECont& c1, const vector<ECont>& c2, TRet&& ret){
	ECont ac = c1;
	for (int i=0; i<c2.size(); ++i){
		ac.insert(ac.end(), c2[i].begin(), c2[i].end());
	}
	ECont ae = ret.alledges();
	HM2D::ECol::Algos::AssignBTypes(ac, ae);
	return ret;
}
};

//#define USE_LIBCLIPPER_FOR_CLIPPING

#ifdef USE_LIBCLIPPER_FOR_CLIPPING //using libclipper

//two contours. direction is not taken into account
TRet ci::Intersection(const ECont& c1, const ECont& c2){
	vector<Impl::ClipperPath> p1{c1}, p2{c2};
	return assign_btypes(
		c1, c2,
		Impl::ClipperPath::Intersect(p1, p2, false, false));
}

TRet ci::Union(const ECont& c1, const ECont& c2){
	vector<Impl::ClipperPath> p1{c1}, p2{c2};
	return assign_btypes(
		c1, c2,
		Impl::ClipperPath::Union(p1, p2, false, false));
}

TRet ci::Difference(const ECont& c1, const ECont& c2){
	vector<Impl::ClipperPath> p1{c1}, p2{c2};
	return assign_btypes(
		c1, c2,
		Impl::ClipperPath::Substruct(p1, p2, false, false));
}


//tree and contour: contours direction is not taken into account
TRet ci::Intersection(const ETree& c1, const ECont& c2){
	_THROW_NOT_IMP_;
}

TRet ci::Union(const ETree& c1, const ECont& c2){
	vector<Impl::ClipperPath> p1, p2{c2};
	for (auto& n: c1.nodes){
		p1.push_back(Impl::ClipperPath(n->contour));
	}
	return assign_btypes(
		c1.alledges(), c2,
		Impl::ClipperPath::Union(p1, p2, true, false));
}

TRet ci::Difference(const ETree& c1, const ECont& c2){
	vector<Impl::ClipperPath> p1, p2{c2};
	for (auto& n: c1.nodes){
		p1.push_back(Impl::ClipperPath(n->contour));
	}
	return assign_btypes(
		c1.alledges(), c2,
		Impl::ClipperPath::Substruct(p1, p2, true, false));
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
	return assign_btypes(
		cont,
		Impl::ClipperPath::Union(p1, zero, false, false));
}

TRet ci::Difference(const ETree& c1, const vector<ECont>& cont){
	if (cont.size() == 0) return TRet::DeepCopy(c1);
	vector<Impl::ClipperPath> p1, p2;
	for (auto& c: c1.nodes) p1.push_back(Impl::ClipperPath(c->contour));
	for (auto& c: cont) p2.push_back(Impl::ClipperPath(c));
	return assign_btypes(
		c1.alledges(), cont,
		Impl::ClipperPath::Substruct(p1, p2, true, false));
}

void ci::Heal(TRet& c1){
	//make union to connect adjacent areas if they have appeared
	if (c1.nodes.size() > 1){
		vector<Impl::ClipperPath> p1, zero;
		for (auto& c: c1.nodes) p1.push_back(Impl::ClipperPath(c->contour));
		c1 = Impl::ClipperPath::Union(p1, zero, true, false);
	}
	//remove zero length edges
	for (auto& path: c1.nodes){
		auto& cont = path->contour;
		auto pts = HM2D::Contour::CornerPoints(cont);
		if (cont.size() < 3){ 
			//invalid contour. Will be removed.
			cont.resize(2);
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
		auto c = HM2D::Contour::Assembler::Contour1(pts, true);
		path->contour = c;
	}
	//remove invalid contours
	std::set<HM2D::Contour::Tree::TNode*> not_needed;
	for (int i=0; i<c1.nodes.size(); ++i)
		if (c1.nodes[i]->contour.size()<3) not_needed.insert(c1.nodes[i].get());
	for (auto n: not_needed) c1.remove_contour(n);
}

#else //same using gpc

//two contours. direction is not taken into account
TRet ci::Intersection(const ECont& c1, const ECont& c2){
	Impl::GpcTree p1(c1), p2(c2);
	return assign_btypes(
		c1, c2,
		Impl::GpcTree::Intersect(p1, p2).ToContourTree());
}

TRet ci::Union(const ECont& c1, const ECont& c2){
	Impl::GpcTree p1(c1), p2(c2);
	return assign_btypes(
		c1, c2,
		Impl::GpcTree::Union(p1, p2).ToContourTree());
}

TRet ci::Difference(const ECont& c1, const ECont& c2){
	Impl::GpcTree p1(c1), p2(c2);
	return assign_btypes(
		c1, c2,
		Impl::GpcTree::Substract(p1, p2).ToContourTree());
}

//tree and contour: contours direction is not taken into account
TRet ci::Intersection(const ETree& c1, const ECont& c2){
	Impl::GpcTree p1(c1), p2(c2);
	return assign_btypes(
		c1.alledges(), c2,
		Impl::GpcTree::Intersect(p1, p2).ToContourTree());
}

TRet ci::Union(const ETree& c1, const ECont& c2){
	Impl::GpcTree p1(c1), p2(c2);
	return assign_btypes(
		c1.alledges(), c2,
		Impl::GpcTree::Union(p1, p2).ToContourTree());
}

TRet ci::Difference(const ETree& c1, const ECont& c2){
	Impl::GpcTree p1(c1), p2(c2);
	return assign_btypes(
		c1.alledges(), c2,
		Impl::GpcTree::Substract(p1, p2).ToContourTree());
}

TRet ci::Difference(const ECont& c1, const ETree& c2){
	_THROW_NOT_IMP_;
}

//two trees
TRet ci::Intersection(const ETree& c1, const ETree& c2){
	Impl::GpcTree p1(c1), p2(c2);
	return assign_btypes(
		c1.alledges(), c2.alledges(),
		Impl::GpcTree::Intersect(p1, p2).ToContourTree());
}

TRet ci::Union(const ETree& c1, const ETree& c2){
	Impl::GpcTree p1(c1), p2(c2);
	return assign_btypes(
		c1.alledges(), c2.alledges(),
		Impl::GpcTree::Union(p1, p2).ToContourTree());
}

TRet ci::Difference(const ETree& c1, const ETree& c2){
	Impl::GpcTree p1(c1), p2(c2);
	return assign_btypes(
		c1.alledges(), c2.alledges(),
		Impl::GpcTree::Substract(p1, p2).ToContourTree());
}

TRet ci::XOR(const ETree& c1, const ETree& c2){
	Impl::GpcTree p1(c1), p2(c2);
	return assign_btypes(
		c1.alledges(), c2.alledges(),
		Impl::GpcTree::Xor(p1, p2).ToContourTree());
}


//multiple contours operation.
//Direction of each contour is not taken into account.
TRet ci::Intersection(const vector<ECont>& cont){
	_THROW_NOT_IMP_;
}

TRet ci::Union(const vector<ECont>& cont){
	if (cont.size() == 0) return TRet();
	vector<Impl::GpcTree> p1; p1.reserve(cont.size());
	for (auto& c: cont) p1.push_back(Impl::GpcTree(c));
	for (int i=1; i<cont.size(); ++i){
		p1[0] = Impl::GpcTree::Union(p1[0], p1[i]);
	}
	return assign_btypes(
		cont,
		p1[0].ToContourTree());
}

TRet ci::Difference(const ETree& c1, const vector<ECont>& cont){
	if (cont.size() == 0) return Contour::Tree::DeepCopy(c1);
	Impl::GpcTree p1(c1);
	for (auto& c: cont){
		Impl::GpcTree g(c);
		p1 = Impl::GpcTree::Substract(p1, g);
	}
	return assign_btypes(
		c1.alledges(), cont,
		p1.ToContourTree());
}

void ci::Heal(TRet& c1){
	//remove zero length edges
	for (auto& path: c1.nodes){
		auto pts = Contour::CornerPoints(path->contour);
	
		//remove collinear nodes: 180 and 360 angles)
STARTFOR:
		if (pts.size()>2) for (int i=0; i<pts.size(); ++i){
			auto& p0 = pts[(i==0)?pts.size()-1:i-1];
			auto& p1 = pts[i];
			auto& p2 = pts[(i==pts.size()-1)?0:i+1];
			if (fabs(triarea(*p0, *p1, *p2))<geps2){
				pts.erase(pts.begin() + i);
				goto STARTFOR;
			}
		}

		//write simplified contour
		if (pts.size() != path->contour.size()){
			path->contour = HM2D::Contour::Assembler::Contour1(pts, true);
		}
	}
	//remove invalid contours with less then 3 edges
	std::set<int> not_needed;
	for (int i=0; i<c1.nodes.size(); ++i){
		if (c1.nodes[i]->contour.size()<3) not_needed.insert(i);
	}
	if (not_needed.size() > 0){
		aa::remove_entries(c1.nodes, not_needed);
		//updates topology
		c1.update_topology();
	}
}

#endif
