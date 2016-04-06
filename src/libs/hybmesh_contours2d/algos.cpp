#include "algos.hpp"
#include "constructor.hpp"
#include "clipper_core.hpp"

using namespace HMCont2D;
namespace cns = Algos;


// ===================================== Offset implementation
Container<ContourTree> cns::Offset(const Contour& source, double delta, OffsetTp tp){
	Impl::ClipperPath cp(source);
	if (source.is_closed() && Contour::Area(source) < 0) delta = -delta;
	return cp.Offset(delta, tp);
};

Container<Contour> cns::Offset1(const Contour& source, double delta){
	Container<ContourTree> ans;
	if (source.is_closed()) ans = Offset(source, delta, OffsetTp::RC_CLOSED_POLY);
	else ans = Offset(source, delta, OffsetTp::RC_OPEN_ROUND);
	assert(ans.cont_count() == 1);
	return Container<Contour>::DeepCopy(*ans.get_contour(0));
};

// ===================================== Simplifications
ContourTree cns::Simplified(const ContourTree& t1){
	ContourTree ret;
	//copy all contours to ret with simplified structure
	for (auto& n: t1.nodes){
		auto simpcont = Simplified(*n);
		ret.nodes.push_back(std::make_shared<ContourTree::TreeNode>());
		ret.nodes.back()->data = simpcont.data;
	}

	//fill parent
	for (int i=0; i<t1.nodes.size(); ++i){
		auto oldparent = t1.nodes[i]->parent;
		if (oldparent == 0) ret.nodes[i]->parent=0;
		else{
			int fnd=0;
			while (fnd<t1.nodes.size()){
				if (t1.nodes[fnd].get() == oldparent) break;
				else ++fnd;
			}
			assert(fnd<t1.nodes.size());
			ret.nodes[i]->parent = ret.nodes[fnd].get();
		}
	}
	//fill children
	for (auto& c: ret.nodes){
		if (c->parent != 0) c->parent->children.push_back(c.get());
	}
	//fill ret.data
	ret.ReloadEdges();
	return ret;
}

ExtendedTree cns::Simplified(const ExtendedTree& t1){
	ExtendedTree ret;
	//insert simplified closed contour nodes
	auto ct = Simplified(static_cast<ContourTree>(t1));
	ret.nodes.insert(ret.nodes.end(), ct.nodes.begin(), ct.nodes.end());
	//insert simplified open contour nodes
	for (auto& oc: t1.open_contours){
		ret.open_contours.push_back(std::make_shared<Contour>(Simplified(*oc)));
	}
	ret.ReloadEdges();
	return ret;
}

Contour cns::Simplified(const Contour& cont){
	auto p = cont.corner_points();
	return HMCont2D::Constructor::ContourFromPoints(p, cont.is_closed());
}


// ==================================== Crosses and intersections
bool cns::DoIntersect(const Contour& c1, const Contour& c2){
	auto bbox1 = Contour::BBox(c1);
	auto bbox2 = Contour::BBox(c2);
	if (!bbox1.has_common_points(bbox2)) return false;
	auto c = HMCont2D::Clip::Union(c1, c2);
	if (c.cont_count() > 1) return false;
	else return true;
}

bool cns::DoReallyIntersect(const Contour& c1, const Contour& c2){
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

bool cns::DoIntersect(const ContourTree& t1, const Contour& c2){
	auto bbox1 = ContourTree::BBox(t1);
	auto bbox2 = Contour::BBox(c2);
	if (!bbox1.has_common_points(bbox2)) return false;
	auto c = HMCont2D::Clip::Union(t1, c2);
	if (c.cont_count() > t1.cont_count()) return false;
	else return true;
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
	TCoordSet w1;
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
cns::Cross(const Contour& c1, const Contour& c2){
	auto retv = cross_core(c1, c2, true);
	if (retv.size() == 0) return std::make_tuple(false, Point(0,0), 0.0, 0.0);
	else return retv[0];
}

vector<std::tuple<bool, Point, double, double>>
cns::CrossAll(const Contour& c1, const Contour& c2){
	return cross_core(c1, c2, false);
}

vector<int> cns::SortOutPoints(const Contour& t1, const vector<Point>& pnt){
	ContourTree tree;
	tree.AddContour(t1);
	auto ret = cns::SortOutPoints(tree, pnt);
	//if t1 is inner contour
	if (t1.data[0] != tree.nodes[0]->data[0]){
		for (auto& r: ret){
			if (r == OUTSIDE) r = INSIDE;
			else if (r == INSIDE) r = OUTSIDE;
		}
	}
	return ret;
}
vector<int> cns::SortOutPoints(const ContourTree& t1, const vector<Point>& pnt){
	auto ctree = Impl::ClipperTree::Build(HMCont2D::Algos::Simplified(t1));
	return ctree.SortOutPoints(pnt);
}

// =================================== Smoothing
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

Vect cns::SmoothedDirection(const Contour& c, Point* p, int direction, double len){
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
		auto ce = ECollection::FindClosestEdge(c, *p);
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

namespace{
Point smoothed_direction_step(const HMCont2D::Contour& c, double w0, double w_step){
	double w = w0+w_step;
	//adjust w_step
	if (w > 1.0){
		if (!c.is_closed()) return *c.last();
		while (w>1) w -= 1.0;
	}
	return HMCont2D::Contour::WeightPoint(c, w);
}
}
Vect cns::SmoothedDirection2(const Contour& c, const Point *p, int direction, double len_forward, double len_backward){
	//preliminary simplification
	auto cont = HMCont2D::Algos::Simplified(c);
	//descrease lens to half of contour lengths
	double full_len = cont.length();
	if (len_forward > full_len/2) len_forward = full_len/2;
	if (len_backward > full_len/2) len_backward = full_len/2;
	//find points
	double pw = std::get<1>(cont.coord_at(*p));
	Point p1 = smoothed_direction_step(cont, pw, len_forward/full_len);
	cont.ReallyReverse();
	Point p2 = smoothed_direction_step(cont, 1 - pw, len_backward/full_len);
	
	//return zero if all lengths are 0
	if (p1 == p2 && p1 == *p) return Vect(0, 0);

	Vect ret;
	if (p1 != p2) ret = p1 - p2;
	else{
		ret = p1 - *p;
		vecRotate(ret, -M_PI/2);
	}
	if (direction == -1) ret *= -1;
	vecNormalize(ret);
	return ret;
}

