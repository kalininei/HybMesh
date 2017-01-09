#include "cont_assembler.hpp"
#include "algos.hpp"
#include "debug2d.hpp"
#include <stack>
#include "contabs2d.hpp"
#include "treverter2d.hpp"

using namespace HM2D;
using namespace HM2D::Contour;
namespace cns = HM2D::Contour::Assembler;

namespace {

typedef std::list<shared_ptr<Edge>> AList;
typedef std::map<const Point*, AList> AMap;

void sort_angles(const Point* p1, AList& eds){
	if (eds.size() < 3) return;
	std::map<double, shared_ptr<Edge>>  srt;
	std::function<void(double, shared_ptr<Edge>)> add_2_map;
	add_2_map = [&](double a, shared_ptr<Edge> e){
		auto isthere = srt.emplace(a, e);
		if (!isthere.second) add_2_map(a+geps, e);
	};
	for (auto e: eds){
		Point c = *e->sibling(p1) - *p1;
		if (ISZERO(c.x) && ISZERO(c.y)) add_2_map(0, e);
		else add_2_map(atan2(c.y, c.x), e);
	}
	eds.clear();
	for (auto m: srt) eds.push_back(m.second);
}
int calc_more1(const AMap& mp){
	int ret = 0;
	for (auto& m: mp) if (m.second.size()<2) ++ret;
	return ret;
}
AMap::iterator find_largest_list_entry(AMap& mp){
	auto ret = mp.begin();
	for (auto it = mp.begin(); it != mp.end(); ++it){
		if (it->second.size() > ret->second.size()) ret = it;
	}
	return ret;
}

void remove_with_connections(AMap& mp, AMap::iterator pos){
	//remove entry from siblings
	for (auto it = pos->second.begin(); it != pos->second.end(); ++it){
		shared_ptr<Edge> ed = *it;
		const Vertex* sibl = ed->sibling(pos->first).get();
		mp[sibl].remove(ed);
	}
	//remove this section
	mp.erase(pos);
}

void remove_edge_from_amap(AMap& mp, shared_ptr<Edge> ed){
	for (auto& it: mp) it.second.remove(ed);
}


//this splits contour by repeating points
//returns true if input contour was splitted
//and false otherwise. In the latter case no new members to ret is added.
bool get_non_crossed(const EdgeData& cont, vector<EdgeData>& ret){
	VertexData v_allpts = OrderedPoints(cont);
	std::map<Point*, int> m_allpts;
	for (int i=0; i<cont.size(); ++i){
		auto isthere = m_allpts.emplace(v_allpts[i].get(), i);
		if (!isthere.second){
			EdgeData cont1 = cont;
			EdgeData cont2;
			auto it1 = cont.begin() + isthere.first->second;
			auto it2 = cont.begin() + i;
			auto it11 = cont1.begin() + isthere.first->second;
			auto it12 = cont1.begin() + i;
			cont1.erase(it11, it12);
			cont2.insert(cont2.end(), it1, it2);
			ret.push_back(cont2);
			if (!get_non_crossed(cont1, ret)) ret.push_back(cont1);
			return true;
		}
	}
	return false;
}

EdgeData _assemble(AMap& mp, AMap::iterator start){
	EdgeData ret;
	auto take_next = [&](AMap::iterator& it, shared_ptr<Edge>& e)->void{
		const Point* p1 = it->first;
		const Point* p2 = e->sibling(p1).get();
		it = mp.find(p2);
		if (it == mp.end()) return;
		if (it == start || it->second.size() == 1) { it = mp.end(); return; }
		if (it->second.size() == 2){
			if (*it->second.begin() == e) e = *it->second.rbegin();
			else e = *it->second.begin();
			return;
		}
		auto fnd_in_list = std::find(it->second.begin(), it->second.end(), e);
		++fnd_in_list;
		if (fnd_in_list == it->second.end()) fnd_in_list = it->second.begin();
		e = *fnd_in_list;
	};

	AMap::iterator it = start;
	shared_ptr<Edge> ecur = *start->second.begin();
	while (it != mp.end()){
		shared_ptr<Edge> eprev = ecur;
		take_next(it, ecur);
		ret.push_back(eprev);
		remove_edge_from_amap(mp, eprev);
	}
	//go backwards if start point had two connected edges
	//and contour was not closed
	if (IsOpen(ret) && start->second.size() == 1){
		Reverse(ret);
		ecur = *start->second.begin();
		it = start;
		while (it != mp.end()){
			shared_ptr<Edge> eprev = ecur;
			take_next(it, ecur);
			ret.push_back(eprev);
			remove_edge_from_amap(mp, eprev);
		}
		Reverse(ret);
	}
	return ret;
}

void remove_zero_list_entries(AMap& mp){
	auto it = mp.begin();
	while (it != mp.end()){
		if (it->second.size() == 0) it = mp.erase(it);
		else ++it;
	}
}
AMap::iterator find_size1_list_entry(AMap& mp){
	for (auto it=mp.begin(); it!=mp.end(); ++it){
		if (it->second.size() == 1) return it;
	}
	return mp.end();
}

void refill_amap(AMap& mp){
	//copy all edges to set
	std::set<shared_ptr<Edge>> s_edges;
	for (auto m: mp) for (auto e: m.second) s_edges.insert(e);
	//build point entries
	mp.clear();
	for (auto e: s_edges){
		auto f1 = mp.emplace(e->first().get(), AList{e});
		auto f2 = mp.emplace(e->last().get(), AList{e});
		if (!f1.second)  f1.first->second.push_back(e);
		if (!f2.second)  f2.first->second.push_back(e);
	}
	
	//sort angles if needed
	for (auto& m: mp) sort_angles(m.first, m.second);
}

vector<EdgeData> sort_all_contours(vector<EdgeData>& ret, std::map<Edge*, int>& edind){
	for (auto& c: ret){
		if (IsClosed(c)){
			//close contours always outer and start with edge with minimum index
			if (Contour::Area(c) < 0) Reverse(c);
			int mini = edind[c[0].get()]; int mini_ind=0;
			for (int i=1; i<c.size(); ++i){
				int ind = edind[c[i].get()];
				if (ind < mini){
					mini = ind;
					mini_ind = i;
				}
			}
			std::rotate(c.begin(), c.begin() + mini_ind, c.end());
		} else {
			//open contours first edge has lower index than last edge
			if (edind[c[0].get()] > edind[c.back().get()]) Reverse(c);
		}
	}
	//sort contours within ret according to first edge global index
	vector<EdgeData*> v;
	for (auto& c: ret) v.push_back(&c);
	std::sort(v.begin(), v.end(), [&edind](EdgeData* c1, EdgeData* c2){
			return edind[(*c1)[0].get()] < edind[(*c2)[0].get()];
		});
	vector<EdgeData> ret1;
	for (auto cp: v) ret1.push_back(*cp);
	return ret1;
}
}

std::vector<EdgeData> cns::AllContours(const EdgeData& input){
	//1) build point->edges connectivity
	AMap point_edges;
	point_edges[0] = AList(input.begin(), input.end());
	refill_amap(point_edges);

	//3) assembling closed contours starting from points with largest number of connections
	vector<EdgeData> conts;
	AMap open_parts;
	while (point_edges.size() > 0){
		//get rid of single and zero connected point entries
		if (calc_more1(point_edges) > 0) while(1){
			//zero lists
			remove_zero_list_entries(point_edges);
			//list of size = 1
			auto it1 = find_size1_list_entry(point_edges);
			if (it1 == point_edges.end()) break;
			//paste this section to open parts
			open_parts.insert(*it1);
			//remove it from point_edges
			remove_with_connections(point_edges, it1);
		}
		//find largest
		auto cur = find_largest_list_entry(point_edges);
		if (cur == point_edges.end()) break;
		//assemble
		conts.push_back(_assemble(point_edges, cur));
	}

	//4) sort open_parts according to angle
	refill_amap(open_parts);
	while (open_parts.size() > 0){
		auto cur = find_size1_list_entry(open_parts);
		//start assembling from single adjecent vertices
		conts.push_back(_assemble(open_parts, cur));
		remove_zero_list_entries(open_parts);
	}

	//5) split by repeating points
	std::vector<EdgeData> ret;
	for (auto& c: conts){
		if (!get_non_crossed(c, ret)) ret.push_back(c);
	}

	//6) sort contours in order to always obtain the same result with the same given
	//   ecollection edges order
	std::map<Edge*, int> edind;
	int ei = 0;
	for (auto e: input) edind[e.get()] = ei++;
	ret = sort_all_contours(ret, edind);

	return ret;
}

namespace{
struct _DoubleNodePoints{
	typedef std::tuple<int, int, shared_ptr<Vertex>> TBup;
	EdgeData* inp;
	VertexData pheap;
	vector<TBup> tempp;

	_DoubleNodePoints(const EdgeData& input){
		inp = const_cast<EdgeData*>(&input);
		auto pe = HM2D::Connectivity::VertexEdge(*inp);
		for (int i=0; i<pe.size(); ++i) if (pe[i].size() > 2){
			auto p = pe[i].v;
			for (int j=0; j<pe[i].size(); ++j){
				pheap.emplace_back(new Vertex(*p));
				auto edge = (*inp)[pe[i].eind[j]];
				if (edge->first() == p){
					edge->vertices[0] = pheap.back();
					tempp.emplace_back(pe[i].eind[j], 0, p);
				} else {
					assert(edge->last() == p);
					edge->vertices[1] = pheap.back();
					tempp.emplace_back(pe[i].eind[j], 1, p);
				}
			}
		}
	}
	~_DoubleNodePoints(){
		for (auto& it: tempp){
			auto edge = (*inp)[std::get<0>(it)];
			if (std::get<1>(it) == 0)
				edge->vertices[0] = std::get<2>(it);
			else
				edge->vertices[1] = std::get<2>(it);
		}
	}
};
}
vector<EdgeData> cns::SimpleContours(const EdgeData& input){
	//temporary double all multiply (> 2) connected points
	//changes will be reversed at 'tmp' destruction
	auto tmp = _DoubleNodePoints(input);
	//build contours
	vector<EdgeData> ret = AllContours(input);
	//return
	return ret;
}

namespace{
bool contains_point(const EdgeData& ed, const Point* p){
	for (auto e: ed){
		if (e->first().get() == p) return true;
		if (e->last().get() == p) return true;
	}
	return false;
}
EdgeData get_contour_by_pts(const EdgeData& col, const Point* pnt1, const Point* pnt2=0){
	AMap point_edges;
	point_edges[0] = AList(col.begin(), col.end());
	refill_amap(point_edges);
	auto it = point_edges.find(pnt1);
	assert(it != point_edges.end());
	EdgeData ret = _assemble(point_edges, it);
	std::vector<EdgeData> retvec;
	if (get_non_crossed(ret, retvec)){
		std::function<bool(EdgeData&)> fndfunc;
		if (pnt2 == 0) fndfunc = [&](EdgeData& c){ return contains_point(c, pnt1); };
		else fndfunc = [&](EdgeData& c){
			return contains_point(c, pnt1) && contains_point(c, pnt2);
		};
		auto it = std::find_if(retvec.begin(), retvec.end(), fndfunc);
		assert(it != retvec.end());
		ret = *it;
	}
	assert(contains_point(ret, pnt1));
	assert(pnt2 == 0 || contains_point(ret, pnt2));
	//force same result for
	if (IsOpen(ret)){
		if (ret.size() > 1 &&
		    Point::meas(*First(ret), *pnt1) > Point::meas(*Last(ret), *pnt1)){
			Reverse(ret);
		}
	} else {
		if (Area(ret) < 0) Reverse(ret);
	}

	return ret;
}

EdgeData assemble_core(const EdgeData& con, int estart, int eend){
	EdgeData ret;
	if (eend<=con.size()){
		std::copy(con.begin() + estart, con.begin() + eend,
			std::back_inserter(ret));
	} else {
		eend -= con.size();
		std::copy(con.begin() + estart, con.end(), std::back_inserter(ret));
		std::copy(con.begin(), con.begin() + eend, std::back_inserter(ret));
	}
	return ret;
}

}

//Assemble single contour from shattered edges starting from given points of collection edges
EdgeData cns::Contour1(const EdgeData& col, const Point* pnt_start, const Point* pnt_end){
	auto ret = get_contour_by_pts(col, pnt_start, pnt_end);
	if (pnt_end == 0) return ret;
	else return ShrinkContour(ret, pnt_start, pnt_end);
}

//Assemble from another contour
EdgeData cns::ShrinkContour(const EdgeData& con, const Point* pnt_start, const Point* pnt_end){
	assert(IsContour(con));
	//find indicies of pnt_start/pnt_end
	auto op = OrderedPoints(con);
	int i0 = std::find_if(op.begin(), op.end(), [&pnt_start](shared_ptr<Vertex> x)
			{return x.get()==pnt_start;}) - op.begin();
	int i1 = std::find_if(op.begin(), op.end(), [&pnt_end](shared_ptr<Vertex> x)
			{return x.get()==pnt_end;}) - op.begin();
	assert(i0 < op.size() && i1 < op.size());

	if (IsClosed(con)){
		if (i0 == i1) i1 = i0 + con.size();
		if (i1 == 0) i1 = con.size();
	}

	if (i1 >= i0) return assemble_core(con, i0, i1);
	else {
		if (IsClosed(con)) return assemble_core(con, i0, i1 + con.size());
		else{
			//if reversed non-closed
			EdgeData r = assemble_core(con, i1, i0);
			Reverse(r);
			return r;
		}
	}
}

EdgeData cns::Contour1(const VertexData& pnt, bool force_closed){
	EdgeData ret;
	if (pnt.size()<2) return ret;
	for (int i=0; i<(int)pnt.size() - 1; ++i){
		ret.emplace_back(new Edge(pnt[i], pnt[i+1]));
	}
	if (force_closed && pnt[0]!=pnt.back()){
		ret.emplace_back(new Edge(pnt.back(), pnt[0]));
	}
	return ret;
}

vector<EdgeData> cns::QuickSeparate(const EdgeData& ecol){
	auto ee = Connectivity::EdgeEdge(ecol);
	vector<bool> used(ecol.size(), false);
	auto use = [&used, &ee](int ind, vector<int>& lst, std::stack<int>& s){
		if (used[ind] != true){
			used[ind] = true;
			lst.push_back(ind);
			for (int eadj: ee[ind])
				if (used[eadj] == false) s.push(eadj);
		}
	};

	vector<vector<int>> edgeslist;
	while (1){
		int ifnd = std::find(used.begin(), used.end(), false) - used.begin();
		if (ifnd == used.size()) break;
		edgeslist.emplace_back();
		auto& lst =edgeslist.back();
		std::stack<int> su; su.push(ifnd);
		while (su.size() > 0){
			int t = su.top(); su.pop();
			use(t, lst, su);
		}
	}
	if (edgeslist.size() == 1) return {ecol};

	vector<EdgeData> ret(edgeslist.size());
	for (int i=0; i<edgeslist.size(); ++i){
		ret[i].resize(edgeslist[i].size());
		for (int j=0; j<edgeslist[i].size(); ++j){
			ret[i][j] = ecol[edgeslist[i][j]];
		}
	}
	return ret;
}
