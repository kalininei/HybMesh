#include "algos.hpp"
#include "constructor.hpp"
#include "clipper_core.hpp"
#include <unordered_map>

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

ECollection cns::Simplified(const ECollection& ecol, double degree_angle, bool id_nobreak){
	double maxcos = cos(degree_angle/180.*M_PI);
	auto domerge = [&maxcos](Point* p0, Point* p1, Point* p2)->bool{
		double cos = vecDot(*p1-*p0, *p2-*p1) / vecLen(*p1-*p0) / vecLen(*p2-*p1);
		return cos>=maxcos;
	};
	ECollection ret;
	vector<Contour> sc = Assembler::SimpleContours(ecol);
	for (auto& c: sc){
		auto op = c.ordered_points();
		vector<Point*> newp(1, op[0]);
		Point* p0 = op[0];
		for (int i=1; i<op.size()-1; ++i){
			if ( (id_nobreak && c.value(i-1).id != c.value(i).id) ||
			     !domerge(p0, op[i], op[i+1]) ){
				newp.push_back(p0 = op[i]);
			}
		}
		newp.push_back(op.back());
		ret.Unite(Constructor::ContourFromPoints(op));
	}
	return ret;
}

namespace{
struct _TEdgeCrossAnalyser{
	vector<TCoordSet> ecross_set;
	_TEdgeCrossAnalyser(int nedges): ecross_set(nedges) {}
	void add_crosses(Edge* e1, Edge* e2, int i1, int i2){
		double ksieta[2];
		SectCross(*e1->pstart, *e1->pend, *e2->pstart, *e2->pend, ksieta);
		if (ksieta[0] == gbig && ksieta[1] == gbig){
			//if edges are parrallel
			isOnSection(*e1->pstart, *e2->pstart, *e2->pend, ksieta[0]);
			if (ISIN_SS(ksieta[0], 0, 1)) ecross_set[i2].insert(ksieta[0]);
			isOnSection(*e1->pend, *e2->pstart, *e2->pend, ksieta[0]);
			if (ISIN_SS(ksieta[0], 0, 1)) ecross_set[i2].insert(ksieta[0]);
			isOnSection(*e2->pstart, *e1->pstart, *e1->pend, ksieta[0]);
			if (ISIN_SS(ksieta[0], 0, 1)) ecross_set[i1].insert(ksieta[0]);
			isOnSection(*e2->pend, *e1->pstart, *e1->pend, ksieta[0]);
			if (ISIN_SS(ksieta[0], 0, 1)) ecross_set[i1].insert(ksieta[0]);
		} else if (!Edge::AreConnected(*e1, *e2)){
			if (ISEQGREATER(ksieta[0], 0) && ISEQLOWER(ksieta[0], 1) &&
			    ISEQGREATER(ksieta[1], 0) && ISEQLOWER(ksieta[1], 1)){
				//if not connected edges cross
				ecross_set[i1].insert(ksieta[0]);
				ecross_set[i2].insert(ksieta[1]);
			}
		}
	}
	void divide_edges(ShpVector<Edge>& ecol, ShpVector<Point>& pcol){
		for (int i=0; i<ecross_set.size(); ++i) if (ecross_set[i].size()>0){
			auto& st = ecross_set[i];
			auto& edge = *ecol[i];
			if (ISEQ(*st.begin(), 0)) st.erase(st.begin());
			if (st.size() == 0) continue;
			if (ISEQ(*st.rbegin(), 1)) st.erase(std::prev(st.end()));
			if (st.size() == 0) continue;
			Point* p1 = edge.pstart;
			Point* p2 = edge.pend;
			vector<Point*> pa(1, p1);
			for (auto ksi: st){
				pcol.emplace_back(new Point(Point::Weigh(*p1, *p2, ksi)));
				pa.push_back(pcol.back().get());
			}
			pa.push_back(p2);
			for (int i=0; i<pa.size()-1; ++i){
				ecol.emplace_back(new Edge(edge));
				ecol.back()->pstart = pa[i];
				ecol.back()->pend = pa[i+1];
				//edge with equal bounds will be removed in 
				//MergePoints procedure
				edge.pstart = edge.pend = pa[0];
			}
		}
	}
};
}

Container<ECollection> cns::NoCrosses(const ECollection& ecol){
	Container<ECollection> ret;
	Container<ECollection>::DeepCopy(ecol, ret);
	//Find crosses
	BoundingBox area = ECollection::BBox(ret, geps);
	BoundingBoxFinder finder(area, area.maxlen()/20);
	for (auto e: ret) finder.addentry(e->bbox());
	_TEdgeCrossAnalyser ec(ret.size());
	for (int i=0; i<ecol.size(); ++i){
		auto s = finder.suspects(ret.value(i).bbox());
		for (auto ei: s) if (ei!=i){
			ec.add_crosses(ret.pvalue(i), ret.pvalue(ei), i, ei);
		}
	}
	//Part edges
	ec.divide_edges(ret.data, ret.pdata.data);
	//Merge points
	MergePoints(ret);
	DeleteUnusedPoints(ret);
	return ret;
}

// =================================== merge
namespace{
struct _TEdgeSet{
	std::set<std::pair<Point*, Point*>> data;
	vector<bool> added;
	void add_edge(Edge* ed){
		added.resize(added.size() + 1);
		Point* p1 = ed->pstart;
		Point* p2 = ed->pend;
		added.back() = false;
		if (p1 != p2){
			if (p1>p2) std::swap(p1, p2);
			auto er = data.emplace(p1, p2);
			added.back() = er.second;
		}
	}
	bool was_used(int i){ return added[i]; }
};
}
void cns::MergePoints(ECollection& ecol){
	auto ap = ecol.all_points();
	auto pe = ecol.tab_points_edges();
	vector<std::vector<int>> shadows(ap.size());
	vector<bool> isactive(ap.size(), true);
	for (int i=0; i<ap.size(); ++i)
	for (int j=0; j<i; ++j) if (isactive[j]){
		if (Point::meas(*ap[i], *ap[j]) < geps*geps){
			shadows[j].push_back(i);
			isactive[i] = false;
			break;
		}
	}
	for (int i=0; i<ap.size(); ++i)
	for (auto psh: shadows[i]){
		for (auto e: pe[psh]){
			auto edge = ecol.pvalue(e);
			if (edge->pstart == ap[psh])
				edge->pstart = ap[i];
			else
				edge->pend = ap[i];
		}
	}
	_TEdgeSet es;
	for (auto e: ecol) es.add_edge(e.get());
	ShpVector<Edge> newedge;
	for (int i=0; i<ecol.size(); ++i) if (es.was_used(i))
		newedge.push_back(ecol.data[i]);
	std::swap(newedge, ecol.data);
}

void cns::DeleteUnusedPoints(Container<ECollection>& econt){
	auto ap1 = econt.all_points();
	auto& ap2 = econt.pdata.data;
	std::unordered_set<Point*> ap1set(ap1.begin(), ap1.end());
	vector<int> used_points;
	for (int i=0; i<ap2.size(); ++i){
		if (ap1set.find(ap2[i].get()) != ap1set.end()) used_points.push_back(i);
	}
	if (used_points.size() == ap2.size()) return;
	ShpVector<Point> newap2;
	for (int i: used_points) newap2.push_back(ap2[i]);
	std::swap(econt.pdata.data, newap2);
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
			SectCross(*op1[i], *op1[i+1], *op2[j], *op2[j+1], ksieta);
			if (ksieta[0]>-geps && ksieta[0]<1+geps && ksieta[1]>-geps && ksieta[1]<1+geps){
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
	//decrease lens to half of contour lengths
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

