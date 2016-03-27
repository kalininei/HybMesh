#include "cont_assembler.hpp"

using namespace HMCont2D;
namespace cns = Assembler;

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
		const Point* sibl = ed->sibling(pos->first);
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
bool get_non_crossed(const HMCont2D::Contour& cont, vector<HMCont2D::Contour>& ret){
	std::vector<Point*> v_allpts = cont.ordered_points();
	std::map<Point*, int> m_allpts;
	for (int i=0; i<cont.size(); ++i){
		auto isthere = m_allpts.emplace(v_allpts[i], i);
		if (!isthere.second){
			HMCont2D::Contour cont1 = cont;
			HMCont2D::Contour cont2;
			auto it1 = cont.data.begin() + isthere.first->second;
			auto it2 = cont.data.begin() + i;
			auto it11 = cont1.data.begin() + isthere.first->second;
			auto it12 = cont1.data.begin() + i;
			cont1.data.erase(it11, it12);
			cont2.data.insert(cont2.data.end(), it1, it2);
			ret.push_back(cont2);
			if (!get_non_crossed(cont1, ret)) ret.push_back(cont1);
			return true;
		}
	}
	return false;
}

HMCont2D::Contour _assemble(AMap& mp, AMap::iterator start){
	HMCont2D::Contour ret;
	auto take_next = [&](AMap::iterator& it, shared_ptr<Edge>& e)->void{
		const Point* p1 = it->first;
		const Point* p2 = e->sibling(p1);
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
		ret.add_value(eprev);
		remove_edge_from_amap(mp, eprev);
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
		auto f1 = mp.emplace(e->pstart, AList{e});
		auto f2 = mp.emplace(e->pend, AList{e});
		if (!f1.second)  f1.first->second.push_back(e);
		if (!f2.second)  f2.first->second.push_back(e);
	}
	
	//sort angles if needed
	for (auto& m: mp) sort_angles(m.first, m.second);
}

vector<Contour> sort_all_contours(vector<Contour>& ret, std::map<Edge*, int>& edind){
	for (auto& c: ret){
		if (c.is_closed()){
			//close contours always outer and start with edge with minimum index
			if (HMCont2D::Contour::Area(c) < 0) c.Reverse();
			int mini = edind[c.data[0].get()]; int mini_ind=0;
			for (int i=1; i<c.size(); ++i){
				int ind = edind[c.data[i].get()];
				if (ind < mini){
					mini = ind;
					mini_ind = i;
				}
			}
			std::rotate(c.data.begin(), c.data.begin() + mini_ind, c.data.end());
		} else {
			//open contours first edge has lower index than last edge
			if (edind[c.data[0].get()] > edind[c.data.back().get()]) c.Reverse();
		}
	}
	//sort contours within ret according to first edge global index
	vector<Contour*> v;
	for (auto& c: ret) v.push_back(&c);
	std::sort(v.begin(), v.end(), [&edind](Contour* c1, Contour* c2){
			return edind[c1->data[0].get()] < edind[c2->data[0].get()];
		});
	vector<Contour> ret1;
	for (auto cp: v) ret1.push_back(*cp);
	return ret1;
}
}

std::vector<Contour> cns::AllContours(const HMCont2D::ECollection& input){
	//1) build point->edges connectivity
	AMap point_edges;
	point_edges[0] = AList(input.begin(), input.end());
	refill_amap(point_edges);

	//3) assembling closed contours starting from points with largest number of connections
	vector<Contour> conts;
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
	std::vector<Contour> ret;
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


HMCont2D::ExtendedTree cns::ETree(const HMCont2D::ECollection& input){
	HMCont2D::ExtendedTree ret;
	//1) assemble all possible contours
	std::vector<Contour> conts = AllContours(input);
	//2) add them to ret
	std::for_each(conts.begin(), conts.end(), [&ret](Contour& a){ ret.AddContour(a); });
	return ret;
}

namespace{
HMCont2D::Contour get_contour_by_pts(const ECollection& col, const Point* pnt1, const Point* pnt2=0){
	AMap point_edges;
	point_edges[0] = AList(col.begin(), col.end());
	refill_amap(point_edges);
	auto it = point_edges.find(pnt1);
	assert(it != point_edges.end());
	HMCont2D::Contour ret = _assemble(point_edges, it);
	std::vector<HMCont2D::Contour> retvec;
	if (get_non_crossed(ret, retvec)){
		std::function<bool(HMCont2D::Contour&)> fndfunc;
		if (pnt2 == 0) fndfunc = [&](HMCont2D::Contour& c){ return c.contains_point(pnt1); };
		else fndfunc = [&](HMCont2D::Contour& c){
			return c.contains_point(pnt1) && c.contains_point(pnt2);
		};
		auto it = std::find_if(retvec.begin(), retvec.end(), fndfunc);
		assert(it != retvec.end());
		ret = *it;
	}
	assert(ret.contains_point(pnt1));
	assert(pnt2 ==0 || ret.contains_point(pnt2));
	return ret;
}

HMCont2D::Contour assemble_core(const Contour& con, int estart, int eend){
	Contour ret;
	if (eend<=con.size()){
		std::copy(con.begin() + estart, con.begin() + eend,
			std::back_inserter(ret.data));
	} else {
		eend -= con.size();
		std::copy(con.begin() + estart, con.end(), std::back_inserter(ret.data));
		std::copy(con.begin(), con.begin() + eend, std::back_inserter(ret.data));
	}
	return ret;
}

}

HMCont2D::Contour cns::Contour1(const ECollection& col, const Point* pnt_start){
	return get_contour_by_pts(col, pnt_start);
}

//Assemble single contour from shattered edges starting from given points of collection edges
HMCont2D::Contour cns::Contour1(const ECollection& col, const Point* pnt_start, const Point* pnt_end){
	auto ret = get_contour_by_pts(col, pnt_start, pnt_end);
	return Contour1(ret, pnt_start, pnt_end);
}

//Assemble from another contour
HMCont2D::Contour cns::Contour1(const Contour& con, const Point* pnt_start, const Point* pnt_end){
	//find indicies of pnt_start/pnt_end
	auto op = con.ordered_points();
	int i0 = std::find(op.begin(), op.end(), pnt_start) - op.begin();
	int i1 = std::find(op.begin(), op.end(), pnt_end) - op.begin();
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

//assemles for pnt_start in the direction (+-1) til the length of contour
//resulting contour will be longer or equal to givenn len
HMCont2D::Contour cns::Contour1(const Contour& col, const Point* pnt_start, int direction, double len){
	assert(col.contains_point(pnt_start));
	if (col.size() == 1) return Contour(col);
	if (direction == -1){
		Contour col2(col);
		col2.Reverse();
		return Contour1(col2, pnt_start, 1, len);
	}

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

	return Contour1(col, pnt_start, pnt_end);
}
