#include "modcont.hpp"
#include "clipper_core.hpp"
#include "contabs2d.hpp"
#include "assemble2d.hpp"
#include "treverter2d.hpp"
#include "finder2d.hpp"

using namespace HM2D;
using Contour::Tree;
namespace cal = HM2D::Contour::Algos;
namespace eal = HM2D::ECol::Algos;

// ===================================== Simplifications
Tree cal::Simplified(const Tree& t1){
	Tree ret = Tree::DeepCopy(t1);
	//copy all contours to ret with simplified structure
	for (auto& n: ret.nodes){
		n->contour = Simplified(n->contour);
	}

	return ret;
}

EdgeData cal::Simplified(const EdgeData& cont){
	auto p = CornerPoints1(cont);
	EdgeData ret = Assembler::Contour1(p, IsClosed(cont));
	auto op = OrderedPoints1(cont);
	aa::enumerate_ids_pvec(op);
	for (auto& e: ret) e->boundary_type = cont[e->pfirst()->id]->boundary_type;
	return ret;
}

namespace{

vector<int> break_by_angle(const VertexData& points, int istart, int iend,
		const vector<double>& angles, double angle0){
	if (angle0>=180) return {istart, iend};
	vector<double> aplus(iend-istart+1, 0.);
	for (int i=istart+1; i<iend; ++i){
		aplus[i-istart] = (fabs(angles[i]-180.0) + aplus[i-istart-1]);
	}
	aplus.back() = aplus.end()[-2];
	double backangle = aplus.back();
	if (points[istart] == points[iend]){
		backangle += fabs(angles[iend]-180);
	}
	if (ISEQLOWER(aplus.back(), angle0)) return {istart, iend};
	int anglestep_n = std::ceil(backangle / angle0);
	double angle_step = backangle/anglestep_n;
	vector<int> ret {istart};
	for (int i=istart+1; i<iend; ++i){
		int div1 = aplus[i-istart] / angle_step;
		//int div2 = (aplus[i-istart+1] - geps)/ angle_step;
		int div2 = aplus[i-istart+1]/ angle_step;
		if (div2 > div1) ret.push_back(i);
	}
	ret.push_back(iend);
	return ret;
}

auto to_keep(Edge* e0, Edge* e1, double angle, double a0, bool bt_nobreak,
		const VertexData& keep)->bool{
	if (bt_nobreak && e0->boundary_type != e1->boundary_type) return true;
	if (angle < 180 - a0 || angle > 180 + a0) return true;
	Vertex* v0 = e0->pfirst();
	if (e1->pfirst() != v0 && e1->plast() != v0)
		v0 = e0->plast();
	if (Finder::Contains(keep, v0)) return true;
	return false;
};

void reset_start_for_closed(EdgeData& c, vector<double>& a, double a0, bool bt_nobreak,
		const VertexData& keep){
	if (c.size()<3) return;
	int istart = 0;
	int icur = 0;
	Edge* eprev = c.back().get();
	while (icur != c.size()){
		if (to_keep(eprev, c[icur].get(), a[icur], a0, bt_nobreak, keep)){
			istart = icur;
			break;
		} else {
			eprev = c[icur++].get();
		}
	}
	//revert data
	if (istart != 0){
		std::rotate(c.begin(), c.begin()+istart, c.end());
		a.pop_back();
		std::rotate(a.begin(), a.begin()+istart, a.end());
		a.push_back(a[0]);
	}
};

};

EdgeData eal::Simplified(const EdgeData& ecol, double degree_angle, bool bt_nobreak,
		const VertexData& keep){
	EdgeData ret;
	if (degree_angle < 0) {
		DeepCopy(ecol, ret, 0);
		return ret;
	} else if (degree_angle == 0){ degree_angle = geps/M_PI*180.; }

	vector<EdgeData> sc = Contour::Assembler::SimpleContours(ecol);
	for (auto& c: sc){
		auto op = Contour::OrderedPoints(c);
		vector<double> angles(op.size(), 0);
		for (int i=1; i<op.size()-1; ++i){
			angles[i] = Angle(*op[i-1], *op[i], *op[i+1])/M_PI*180;
		}
		if (Contour::IsClosed(c)){
			angles[0] = angles.back() = Angle(*op.end()[-2], *op[0], *op[1])/M_PI*180;
			reset_start_for_closed(c, angles, degree_angle, bt_nobreak, keep);
			op = Contour::OrderedPoints(c);
		}
		vector<int> significant_points(1, 0);
		//1 break at ids and sharp angles
		for (int i=1; i<op.size()-1; ++i){
			if (to_keep(c[i-1].get(), c[i].get(), angles[i], degree_angle, 
						bt_nobreak, keep))
				significant_points.push_back(i);
		}
		significant_points.push_back(op.size()-1);
		if (significant_points.size() == op.size()){
			DeepCopy(c, ret, 0);
			continue;
		}
		//2 analyse each section between significant points
		for (int i=0; i<significant_points.size()-1; ++i){
			int i0 = significant_points[i];
			int i1 = significant_points[i+1];
			if (i0+1 == i1){
				ret.push_back(c[i0]);
				continue;
			}
			vector<int> bba = break_by_angle(op, i0, i1, angles, degree_angle);
			for (int j=0; j<bba.size()-1; ++j){
				int j1 = bba[j];
				int j2 = bba[j+1];
				ret.push_back(std::make_shared<Edge>(op[j1], op[j2]));
				ret.back()->boundary_type = c[j]->boundary_type;
			}
		}
	}
	return ret;
}

namespace{
struct _TEdgeCrossAnalyser{
	vector<TCoordSet> ecross_set;
	_TEdgeCrossAnalyser(int nedges): ecross_set(nedges) {}
	void add_crosses(Edge* e1, Edge* e2, int i1, int i2){
		double ksieta[2];
		SectCross(*e1->first(), *e1->last(), *e2->first(), *e2->last(), ksieta);
		if (ksieta[0] == gbig && ksieta[1] == gbig){
			//if edges are parrallel
			isOnSection(*e1->first(), *e2->first(), *e2->last(), ksieta[0]);
			if (ISIN_NN(ksieta[0], 0, 1)) ecross_set[i2].insert(ksieta[0]);
			isOnSection(*e1->last(), *e2->first(), *e2->last(), ksieta[0]);
			if (ISIN_NN(ksieta[0], 0, 1)) ecross_set[i2].insert(ksieta[0]);
			isOnSection(*e2->first(), *e1->first(), *e1->last(), ksieta[0]);
			if (ISIN_NN(ksieta[0], 0, 1)) ecross_set[i1].insert(ksieta[0]);
			isOnSection(*e2->last(), *e1->first(), *e1->last(), ksieta[0]);
			if (ISIN_NN(ksieta[0], 0, 1)) ecross_set[i1].insert(ksieta[0]);
		} else if (!e1->connected_to(*e2)){
			if (ISEQGREATER(ksieta[0], 0) && ISEQLOWER(ksieta[0], 1) &&
			    ISEQGREATER(ksieta[1], 0) && ISEQLOWER(ksieta[1], 1)){
				//if not connected edges cross
				ecross_set[i1].insert(ksieta[0]);
				ecross_set[i2].insert(ksieta[1]);
			}
		}
	}
	void divide_edges(ShpVector<Edge>& ecol){
		for (int i=0; i<ecross_set.size(); ++i) if (ecross_set[i].size()>0){
			auto& st = ecross_set[i];
			auto& edge = *ecol[i];
			auto p1 = edge.first();
			auto p2 = edge.last();
			if (ISEQ(*st.begin(), 0)) st.erase(st.begin());
			if (st.size() == 0) continue;
			if (ISEQ(*st.rbegin(), 1)) st.erase(std::prev(st.end()));
			if (st.size() == 0) continue;
			VertexData pa(1, p1);
			for (auto ksi: st){
				pa.push_back(std::make_shared<Vertex>(Point::Weigh(*p1, *p2, ksi)));
			}
			pa.push_back(p2);
			for (int i=0; i<pa.size()-1; ++i){
				ecol.emplace_back(new Edge(edge));
				ecol.back()->vertices[0] = pa[i];
				ecol.back()->vertices[1] = pa[i+1];
			}
			//edge with equal bounds will be removed in 
			//MergePoints procedure
			edge.vertices[0] = edge.vertices[1] = pa[0];
		}
	}
};
}

EdgeData eal::NoCrosses(const EdgeData& ecol){
	EdgeData ret;
	DeepCopy(ecol, ret);
	//Find crosses
	BoundingBox area = BBox(ret, geps);
	BoundingBoxFinder finder(area, area.maxlen()/20);
	for (auto e: ret) finder.addentry(BoundingBox(*e->first(), *e->last()));
	_TEdgeCrossAnalyser ec(ret.size());
	for (int i=0; i<ecol.size(); ++i){
		auto s = finder.suspects(BoundingBox(*ret[i]->first(), *ret[i]->last()));
		for (auto ei: s) if (ei!=i){
			ec.add_crosses(ret[i].get(), ret[ei].get(), i, ei);
		}
	}
	//Part edges
	ec.divide_edges(ret);
	//Merge points
	MergePoints(ret);
	return ret;
}

// =================================== merge
namespace{
struct _TEdgeSet{
	std::set<std::pair<Point*, Point*>> data;
	vector<bool> added;
	void add_edge(Edge* ed){
		added.resize(added.size() + 1);
		Point* p1 = ed->first().get();
		Point* p2 = ed->last().get();
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
void eal::MergePoints(EdgeData& ecol){
	auto pe = Connectivity::VertexEdge(ecol);
	vector<std::vector<int>> shadows(pe.size());
	vector<bool> isactive(pe.size(), true);
	for (int i=0; i<pe.size(); ++i)
	for (int j=0; j<i; ++j) if (isactive[j]){
		if (*pe[i].v == *pe[j].v){
			shadows[j].push_back(i);
			isactive[i] = false;
			break;
		}
	}
	for (int i=0; i<pe.size(); ++i)
	for (auto psh: shadows[i]){
		for (auto ei: pe[psh].eind){
			auto edge = ecol[ei];
			if (edge->first() == pe[psh].v)
				edge->vertices[0] = pe[i].v;
			else
				edge->vertices[1] = pe[i].v;
		}
	}
	_TEdgeSet es;
	for (auto e: ecol) es.add_edge(e.get());
	ShpVector<Edge> newedge;
	for (int i=0; i<ecol.size(); ++i) if (es.was_used(i))
		newedge.push_back(ecol[i]);
	std::swap(newedge, ecol);
}

void cal::RemovePoints(EdgeData& data, vector<int> ipnt){
	bool is_closed = Contour::IsClosed(data);
	std::sort(ipnt.begin(), ipnt.end());
	auto f = std::upper_bound(ipnt.begin(), ipnt.end(), data.size()+1);
	ipnt.resize(f-ipnt.begin());
	if (ipnt.size() == 0) return;
	if (is_closed && ipnt.back() == data.size()) ipnt.insert(ipnt.begin(), 0);
	ipnt.resize(std::unique(ipnt.begin(), ipnt.end()) - ipnt.begin());
	auto connect_edges = [](Edge* ed0, Edge* ed1)->shared_ptr<Edge>{
		shared_ptr<Edge> newe(new Edge(*ed0));
		if (ed0->last() == ed1->first()){
			newe->vertices[1] = ed1->last();
		} else if (ed0->last() == ed1->last()){
			newe->vertices[1] = ed1->first();
		} else if (ed0->first() == ed1->first()){
			newe->vertices[0] = ed1->last();
		} else {
			newe->vertices[0] = ed1->first();
		}
		return newe;
	};
	while (ipnt.size()>0){
		int in = ipnt.back();
		ipnt.resize(ipnt.size()-1);
		if (is_closed && in == 0){
			shared_ptr<Edge> newe = connect_edges(data.back().get(), data[0].get());
			data.back() = newe;
			data.erase(data.begin());
		} else if (!is_closed && in == data.size()){
			data.resize(data.size()-1);
		} else if (!is_closed && in == 0){
			data.erase(data.begin());
		} else {
			assert(in>0 && in<data.size());
			shared_ptr<Edge> newe = connect_edges(data[in-1].get(), data[in].get());
			data[in-1] = newe;
			data.erase(data.begin()+in);
		}
	}
}

void eal::AssignBTypes(const EdgeData& from, EdgeData& to){
	if (to.size() == 0) return;
	if (from.size() == 0){
		for (auto& e: to) e->boundary_type = 0;
		return;
	}
	//if all from edges have same bt
	vector<int> bt(from.size());
	for (int i=0; i<from.size(); ++i) bt[i] = from[i]->boundary_type;
	if (all_of(bt.begin(), bt.end(), [&bt](int a){ return a == bt[0]; })){
		for (auto& e: to) e->boundary_type = bt[0];
		return;
	}
	//using custom Finder::ClosestEdge implementation
	//to use epsilon comparison instead of '<'. 
	//This is done to guarantee that amoung all equal distanced 'from'
	//edges the first one will be taken.
	for (int ei=0; ei<to.size(); ++ei){
		HM2D::Edge& e = *to[ei];
		Point p = e.center();
		double mindist = std::numeric_limits<double>::max();
		int imin = -1;

		for (int i=0; i<from.size(); ++i){
			double meas = Point::meas_section(p, *from[i]->pfirst(), *from[i]->plast());
			if (meas < geps * geps){
				imin = i;
				break;
			}
			meas = sqrt(meas);
			if (meas < mindist - geps){
				imin = i;
				mindist = meas;
			}
		}
		e.boundary_type = from[imin]->boundary_type;
	}
}

void cal::Reverse(EdgeData& ed){ std::reverse(ed.begin(), ed.end()); }

void cal::AddLastPoint(EdgeData& to, std::shared_ptr<Vertex> p){
	auto p0 = Last(to);
	to.push_back(std::make_shared<Edge>(p0, p));
}


void cal::Connect(EdgeData& to, const EdgeData& from){
	auto self0 = First(to), self1 = Last(to);
	auto target0 = First(from), target1 = Last(from);
	//choosing option for unition
	if (to.size() == 0 || from.size() == 0 ) goto COPY12;
	else if (self0 == self1 || target0 == target1) goto THROW;
	else if (from.size() == 1 &&
		(from[0]->first() == self1 || from[0]->last() == self1)) goto COPY12;
	else if (from.size() == 1 &&
		(from[0]->first() == self0 || from[0]->last() == self0)) goto COPY03;
	//try to add new contour to the end of current
	else if (self1 == target0) goto COPY12;
	else if (self1 == target1) goto NEED_REVERSE;
	//if failed try to add before the start
	else if (self0 == target1) goto COPY03;
	else if (self0 == target0) goto NEED_REVERSE;
	else goto THROW;

COPY03:
	{
		to.insert(to.begin(), from.begin(), from.end());
		return;
	}
COPY12:
	{
		to.insert(to.end(), from.begin(), from.end());
		return;
	}
NEED_REVERSE:
	{
		EdgeData tmp(from);
		Reverse(tmp);
		return Connect(to, tmp);
	}
THROW:
	{
		throw std::runtime_error("Impossible to unite non-connected contours");
	}
}

void cal::SplitEdge(EdgeData& cont, int iedge, const vector<Point>& pts){
	if (pts.size() == 0) return;
	assert(IsContour(cont));
	bool iscor = CorrectlyDirectedEdge(cont, iedge);
	if (!iscor) cont[iedge]->reverse();

	VertexData pp(pts.size()+2);
	pp[0] = cont[iedge]->vertices[0];
	pp.back() = cont[iedge]->vertices[1];
	for (int i=1; i<pp.size()-1; ++i){
		pp[i] = std::make_shared<Vertex>(pts[i-1]);
	}

	cont[iedge]->vertices[1] = pp[1];
	EdgeData newed(pts.size());
	for (int i=0; i<pts.size(); ++i){
		newed[i] = std::make_shared<Edge>(*cont[iedge]);
		newed[i]->vertices[0] = pp[i+1];
		newed[i]->vertices[1] = pp[i+2];
	}
	cont.insert(cont.begin()+iedge+1, newed.begin(), newed.end());

	if (!iscor){
		for (auto it=cont.begin()+iedge; it!=cont.begin()+iedge+pts.size()+1; ++it){
			(*it)->reverse();
		}
	}
};

std::tuple<bool, shared_ptr<Vertex>>
cal::GuaranteePoint(EdgeData& ed, const Point& p){
	std::tuple<bool, shared_ptr<Vertex>> ret;

	auto ce = HM2D::Finder::ClosestEdge(ed, p);
	if (std::get<0>(ce)<0){
		std::get<0>(ret) = false;
		return ret;
	}
	Edge* e = ed[std::get<0>(ce)].get();
	double elen = e->length();
	double len2 = std::get<2>(ce)*elen;
	if (ISZERO(len2)){
		std::get<0>(ret) = false;
		std::get<1>(ret) = e->first();
	} else if (ISZERO(elen-len2)){
		std::get<0>(ret) = false;
		std::get<1>(ret) = e->last();
	} else {
		auto pnew = std::make_shared<Vertex>(Point::Weigh(
			*e->first(), *e->last(), std::get<2>(ce)));
		auto e1 = std::make_shared<Edge>(*e);
		auto e2 = std::make_shared<Edge>(*e);
		e1->vertices[1] = pnew;
		e2->vertices[0] = pnew;
		if (CorrectlyDirectedEdge(ed, std::get<0>(ce))){
			ed[std::get<0>(ce)] = e1;
			ed.insert(ed.begin()+std::get<0>(ce)+1, e2);
		} else {
			ed[std::get<0>(ce)] = e2;
			ed.insert(ed.begin()+std::get<0>(ce)+1, e1);
		}
		std::get<0>(ret) = true;
		std::get<1>(ret) = pnew;
	}
	return ret;
}
