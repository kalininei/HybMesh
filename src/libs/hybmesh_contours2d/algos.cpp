#include "algos.hpp"
#include "clipper_core.hpp"
#include "contabs2d.hpp"
#include "cont_assembler.hpp"
#include "treverter2d.hpp"
#include "finder2d.hpp"

using namespace HM2D;
using Contour::Tree;
namespace cal = HM2D::Contour::Algos;
namespace eal = HM2D::ECol::Algos;

// ===================================== Offset implementation
Tree cal::Offset(const EdgeData& source, double delta, OffsetTp tp){
	Impl::ClipperPath cp(source);
	if (IsClosed(source) && Contour::Area(source) < 0) delta = -delta;
	return cp.Offset(delta, tp);
};

EdgeData cal::Offset1(const EdgeData& source, double delta){
	Contour::Tree ans;
	if (IsClosed(source)) ans = Offset(source, delta, OffsetTp::RC_CLOSED_POLY);
	else ans = Offset(source, delta, OffsetTp::RC_OPEN_ROUND);
	assert(ans.nodes.size() == 1);
	return ans.nodes[0]->contour;
};


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
	auto p = CornerPoints(cont);
	return Assembler::Contour1(p, IsClosed(cont));
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
		int div2 = (aplus[i+1-istart] - geps)/ angle_step;
		if (div2 > div1) ret.push_back(i);
	}
	ret.push_back(iend);
	return ret;
}

}

EdgeData eal::Simplified(const EdgeData& ecol, double degree_angle, bool id_nobreak){
	EdgeData ret;
	if (degree_angle < 0) {
		DeepCopy(ecol, ret, 0);
		return ret;
	}

	vector<EdgeData> sc = Contour::Assembler::SimpleContours(ecol);
	for (auto& c: sc){
		auto op = Contour::OrderedPoints(c);
		vector<double> angles(op.size(), 0);
		for (int i=1; i<op.size()-1; ++i){
			angles[i] = Angle(*op[i-1], *op[i], *op[i+1])/M_PI*180;
		}
		if (Contour::IsClosed(c)){
			angles[0] = angles.back() = Angle(*op.end()[-2], *op[0], *op[1])/M_PI*180;
		}
		vector<int> significant_points(1, 0);
		//1 break at ids and sharp angles
		for (int i=1; i<op.size()-1; ++i){
			if ( (id_nobreak && c[i-1]->id != c[i]->id) ||
			      angles[i] < 180 - degree_angle ||
			      angles[i] > 180 + degree_angle)
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
				ret.back()->id = c[j]->id;
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
			if (ISEQ(*st.begin(), 0)) st.erase(st.begin());
			if (st.size() == 0) continue;
			if (ISEQ(*st.rbegin(), 1)) st.erase(std::prev(st.end()));
			if (st.size() == 0) continue;
			auto p1 = edge.first();
			auto p2 = edge.last();
			VertexData pa(1, p1);
			for (auto ksi: st){
				pa.push_back(std::make_shared<Vertex>(Point::Weigh(*p1, *p2, ksi)));
			}
			pa.push_back(p2);
			for (int i=0; i<pa.size()-1; ++i){
				ecol.emplace_back(new Edge(edge));
				ecol.back()->first() = pa[i];
				ecol.back()->last() = pa[i+1];
				//edge with equal bounds will be removed in 
				//MergePoints procedure
				edge.first() = edge.last() = pa[0];
			}
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

vector<int> cal::SortOutPoints(const EdgeData& t1, const vector<Point>& pnt){
	Contour::Tree tree;
	tree.add_contour(t1);
	auto ret = cal::SortOutPoints(tree, pnt);
	//if t1 is inner contour
	if (Contour::Area(t1) < 0){
		for (auto& r: ret){
			if (r == OUTSIDE) r = INSIDE;
			else if (r == INSIDE) r = OUTSIDE;
		}
	}
	return ret;
}
vector<int> cal::SortOutPoints(const Tree& t1, const vector<Point>& pnt){
	std::function<void(const ShpVector<Tree::TNode>&, Point&, int&)>
	lvwithin = [&lvwithin](const ShpVector<Tree::TNode>& lv, Point& p, int& a){
		int indwithin = -1;
		for (int i=0; i<lv.size(); ++i){
			int rs = Contour::Finder::WhereIs(lv[i]->contour, p);
			if (rs == BOUND){ a = -2; return; }
			else if (rs == INSIDE) { indwithin = i; break; }
		}
		if (indwithin == -1) return;
		else {
			ShpVector<Tree::TNode> newlv;
			for (auto w: lv[indwithin]->children) newlv.push_back(w.lock());
			return lvwithin(newlv, p, ++a);
		}
	};

	vector<int> ret;
	auto roots = t1.roots();
	for (auto p: pnt){
		int maxlevel = -1;
		lvwithin(roots, p, maxlevel);
		if (maxlevel == -2) ret.push_back(BOUND);
		else if (maxlevel == -1) ret.push_back(OUTSIDE);
		else if (maxlevel % 2 == 0) ret.push_back(INSIDE);
		else ret.push_back(OUTSIDE);
	}

	return ret;
}

namespace{
Point smoothed_direction_step(const EdgeData& c, double w0, double w_step){
	double w = w0+w_step;
	//adjust w_step
	if (w > 1.0){
		if (Contour::IsOpen(c)) return *Contour::Last(c);
		while (w>1) w -= 1.0;
	}
	return Contour::WeightPoint(c, w);
}
}
Vect cal::SmoothedDirection2(const EdgeData& c, const Point* p, int direction, double len_forward, double len_backward){
	//preliminary simplification
	auto cont = cal::Simplified(c);
	//decrease lens to half of contour lengths
	double full_len = Length(cont);
	if (len_forward > full_len/2) len_forward = full_len/2;
	if (len_backward > full_len/2) len_backward = full_len/2;
	//find points
	double pw = std::get<1>(Contour::CoordAt(cont, *p));
	Point p1 = smoothed_direction_step(cont, pw, len_forward/full_len);
	Contour::R::ReallyRevert::Permanent(cont);
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

void cal::RemovePoints(EdgeData& data, vector<int> ipnt){
	std::sort(ipnt.begin(), ipnt.end());
	bool is_closed = Contour::IsClosed(data);
	ipnt.resize(std::unique(ipnt.begin(), ipnt.end()) - ipnt.begin());
	while (ipnt.size()>0){
		int in = ipnt.back();
		ipnt.resize(ipnt.size()-1);

		if (is_closed && in == 0){
			auto ed0 = data.back();
			auto ed1 = data[0];
			shared_ptr<Edge> newe(new Edge(*ed0));
			newe->vertices[1] = (ed1->first() == ed0->last()) ? ed1->last()
			                                                  : ed1->first();
			data.back() = newe;
			data.erase(data.begin());
		} else if (!is_closed && in == data.size()){
			data.resize(data.size()-1);
		} else if (!is_closed && in == 0){
			data.erase(data.begin());
		} else {
			assert(in>0 && in<data.size());
			auto ed0 = data[in-1];
			auto ed1 = data[in];
			shared_ptr<Edge> newe(new Edge(*ed0));
			newe->vertices[1] = (ed1->first() == ed0->last()) ? ed1->last()
			                                                  : ed1->first();
			data[in-1] = newe;
			data.erase(data.begin()+in);
		}
	}
}
