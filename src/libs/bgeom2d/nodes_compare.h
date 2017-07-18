#ifndef POINTS_COMPARE_H
#define POINTS_COMPARE_H

#include "bgeom2d.h"

struct DWRef_{ public: double ref; };

inline bool dwref_cmp1(const DWRef_& a, const DWRef_& b){ return a.ref<b.ref; }
inline bool dwref_cmp2(const DWRef_& a, double b){ return a.ref<b; }
inline bool dwref_cmp3(double a, const DWRef_& b){ return a<b.ref; }

// ==================== 2D set
template<class TPoint, class TProc>
class Point2Set{
	struct XData: public DWRef_{
		vector<const TPoint*> subs;

		XData(typename vector<const TPoint*>::iterator i0, typename vector<const TPoint*>::iterator i1): subs(i0, i1){
			auto psorty = [](const TPoint* a, const TPoint* b){
				return TProc::y(*a)<TProc::y(*b);
			};
			ref = TProc::x(**i0) + geps;
			std::sort(subs.begin(), subs.end(), psorty);
		}
	
		//search for equal point in current column
		const TPoint* search_for_point(double x, double y) const{
			auto fnd = std::upper_bound(subs.begin(), subs.end(), y-geps,
				[](double a, const TPoint* b){ return a<TProc::y(*b);});
			while (fnd != subs.end() && TProc::y(**fnd)-y < geps){
				if (fabs(TProc::x(**fnd)-x)<geps) return *fnd;
				++fnd;
			}
			return nullptr;
		}
		vector<const TPoint*> search_for_points(double x, double y, int istart=0) const{
			vector<const TPoint*> ret;
			auto fnd = std::upper_bound(subs.begin() + istart, subs.end(), y-geps,
				[](double a, const TPoint* b){ return a<TProc::y(*b);});
			while (fnd != subs.end() && TProc::y(**fnd)-y < geps){
				if (fabs(TProc::x(**fnd)-x)<geps) ret.push_back(*fnd);
				++fnd;
			}
			return ret;
		}
	};
	vector<XData> sdata;
	const vector<TPoint>* pinput;
	int gind(const TPoint* p) const { return p - &(*pinput)[0]; }
public:
	Point2Set(const vector<TPoint>& input): pinput(&input){
		if (input.size() == 0) return;
		auto psortx = [](const TPoint* a, const TPoint* b){
			return TProc::x(*a)<TProc::x(*b);
		};
		vector<const TPoint*> vp(input.size());
		for (int i=0; i<input.size(); ++i) vp[i] = &input[i];
		std::sort(vp.begin(), vp.end(), psortx);
		auto it=vp.begin();
		while (it != vp.end()){
			auto itend = it+1;
			while (itend != vp.end() && TProc::x(**itend) - TProc::x(**it) < 2*geps) ++itend;
			sdata.push_back(XData(it, itend));
			it = itend;
		}
	}
	//-> vector of <secondary, base> points.
	//   All secondary points will be removed from the structure.
	std::vector<std::pair<int, int>> verbose_unique(){
		std::vector<bool> is_unique(pinput->size(), true);
		std::vector<std::pair<int, int>> ret;
		std::vector<bool> rows_to_resize(sdata.size(), false);

		auto add_doubled_point = [&](int from, int to, int irow){
			if (is_unique[from] == false) return;
			is_unique[from] = false;
			ret.push_back(std::make_pair(from, to));
			rows_to_resize[irow] = true;
		};

		//calculate doubled points
		for (int i=0; i<sdata.size(); ++i){
			for (int j=0; j<sdata[i].subs.size(); ++j){
				const TPoint* p = sdata[i].subs[j];
				int pindex = gind(p);
				if (!is_unique[pindex]) continue;
				//current column
				for (auto p2: sdata[i].search_for_points(TProc::x(*p), TProc::y(*p), j+1)){
					add_doubled_point(gind(p2), pindex, i);
				}
				//next column
				if (TProc::x(*p)>sdata[i].ref && sdata[i+1].ref-TProc::x(*p)<2*geps)
				for (auto p2: sdata[i+1].search_for_points(TProc::x(*p), TProc::y(*p), 0)){
					add_doubled_point(gind(p2), pindex, i+1);
				}
			}
		}
		if (ret.size() == 0) return ret;

		//remove doubled points
		for (int i=0; i<sdata.size(); ++i) if (rows_to_resize[i]){
			vector<const TPoint*> subs2; subs2.reserve(sdata[i].subs.size());
			for (int j=0; j<sdata[i].subs.size(); ++j){
				if (is_unique[gind(sdata[i].subs[j])]) subs2.push_back(sdata[i].subs[j]);
			}
			std::swap(subs2, sdata[i].subs);
		}

		return ret;
	}

	int find(double x, double y) const{
		auto xfnd = std::upper_bound(sdata.begin(), sdata.end(), x, dwref_cmp3);
		if (xfnd != sdata.end() && xfnd->ref-x < 2*geps){
			const TPoint* ret = xfnd->search_for_point(x, y);
			if (ret != nullptr) return gind(ret);
		}
		if (xfnd != sdata.begin() && x-(xfnd-1)->ref < 2*geps){
			const TPoint* ret = (xfnd-1)->search_for_point(x, y);
			if (ret != nullptr) return gind(ret);
		}
		return -1;
	}
};

template<class TPoint>
struct PointerPoint2Access{
	static double x(const TPoint& n) {return n->x;}
	static double y(const TPoint& n) {return n->y;}
};

template<class TPoint>
struct DirectPoint2Access{
	static double x(const TPoint& n) {return n.x;}
	static double y(const TPoint& n) {return n.y;}
};

template<class TPoint>
Point2Set<TPoint, PointerPoint2Access<TPoint>> point2_set_pointers(const std::vector<TPoint>& input){
	return Point2Set<TPoint, PointerPoint2Access<TPoint>>(input);
}

template<class TPoint>
Point2Set<TPoint, DirectPoint2Access<TPoint>> point2_set_coords(const std::vector<TPoint>& input){
	return Point2Set<TPoint, DirectPoint2Access<TPoint>>(input);
}

// ==================== 3D set
template<class TPoint, class TProc>
class Point3Set{
	struct YData: public DWRef_{
		vector<const TPoint*> subs;
		YData(typename vector<const TPoint*>::iterator i0, typename vector<const TPoint*>::iterator i1):subs(i0, i1){
			auto psortz = [](const TPoint* a, const TPoint* b){
				return TProc::z(*a)<TProc::z(*b);
			};
			ref = TProc::y(**i0) + geps;
			std::sort(subs.begin(), subs.end(), psortz);
		}
	
		//search for equal point in current column
		const TPoint* search_for_point(double x, double y, double z) const{
			auto fnd = std::upper_bound(subs.begin(), subs.end(), z-geps,
				[](double a, const TPoint* b){ return a<TProc::z(*b);});
			while (fnd != subs.end() && TProc::z(**fnd)-z<geps){
				if (fabs(TProc::x(**fnd)-x)<geps && fabs(TProc::y(**fnd)-y)<geps) return *fnd;
				++fnd;
			}
			return nullptr;
		}
	};
	struct XData: public DWRef_{
		vector<YData> subs;

		XData(typename vector<const TPoint*>::iterator i0, typename vector<const TPoint*>::iterator i1){
			auto psorty = [](const TPoint* a, const TPoint* b){
				return TProc::y(*a)<TProc::y(*b);
			};
			ref = TProc::x(**i0)+geps;
			std::sort(i0, i1, psorty);
			auto it = i0;
			while (it != i1){
				auto itend = it+1;
				while (itend != i1 && TProc::y(**itend)-TProc::y(**it) < 2*geps) ++itend;
				subs.push_back(YData(it, itend));
				it = itend;
			}
		}
		//search for equal point in current column
		const TPoint* search_for_point(double x, double y, double z) const{
			auto fnd = std::upper_bound(subs.begin(), subs.end(), y, dwref_cmp3);
			if (fnd != subs.end() && fnd->ref - y < 2*geps){
				const TPoint* ret = fnd->search_for_point(x, y, z);
				if (ret != nullptr) return ret;
			}
			if (fnd != subs.begin() && y-(fnd-1)->ref < 2*geps){
				const TPoint* ret = (fnd-1)->search_for_point(x, y, z);
				if (ret != nullptr) return ret;
			}
			return nullptr;
		}

	};
	vector<XData> sdata;
	const vector<TPoint>* pinput;
	int gind(const TPoint* p) const { return p-&(*pinput)[0];}
public:
	Point3Set(const std::vector<TPoint>& input): pinput(&input){
		if (input.size() == 0) return;
		auto psortx = [](const TPoint* a, const TPoint* b){
			return TProc::x(*a)<TProc::x(*b);
		};
		//XData fill
		vector<const TPoint*> vp(input.size());
		for (int i=0; i<input.size(); ++i) vp[i] = &input[i];
		std::sort(vp.begin(), vp.end(), psortx);
		auto it=vp.begin();
		while (it != vp.end()){
			auto itend=it+1;
			while (itend != vp.end() && TProc::x(**itend)-TProc::x(**it) < 2*geps) ++itend;
			sdata.push_back(XData(it, itend));
			it = itend;
		}
	}
	//-> vector of <secondary, base> points.
	//   All secondary points will be removed from the structure.
	std::vector<std::pair<int, int>> verbose_unique(){
		_THROW_NOT_IMP_;
	}

	int find(double x, double y, double z) const{
		auto xfnd = std::upper_bound(sdata.begin(), sdata.end(), x, dwref_cmp3);
		if (xfnd != sdata.end() && xfnd->ref-x < 2*geps){
			const TPoint* ret = xfnd->search_for_point(x, y, z);
			if (ret != nullptr) return gind(ret);
		}
		if (xfnd != sdata.begin() && x-(xfnd-1)->ref < 2*geps){
			const TPoint* ret = (xfnd-1)->search_for_point(x, y, z);
			if (ret != nullptr) return gind(ret);
		}
		return -1;
	}
};

template<class TPoint>
struct PointerPoint3Access{
	static double x(const TPoint& n) {return n->x;}
	static double y(const TPoint& n) {return n->y;}
	static double z(const TPoint& n) {return n->z;}
};

template<class TPoint>
struct DirectPoint3Access{
	static double x(const TPoint& n) {return n.x;}
	static double y(const TPoint& n) {return n.y;}
	static double z(const TPoint& n) {return n.z;}
};

template<class TPoint>
Point3Set<TPoint, PointerPoint3Access<TPoint>> point3_set_pointers(const std::vector<TPoint>& input){
	return Point3Set<TPoint, PointerPoint3Access<TPoint>>(input);
}

template<class TPoint>
Point3Set<TPoint, DirectPoint3Access<TPoint>> point3_set_coords(const std::vector<TPoint>& input){
	return Point3Set<TPoint, DirectPoint3Access<TPoint>>(input);
}

// ==================== 2D point map
template<class D>
class CoordinateMap2D{
	struct XData{
		double realx, realy;
		D data;
	};
	//1 - reference x coordinate (lies in geps neighborhood of real x)
	//2 - real x, y coordinates and data
	//---------1-------------------2
	std::map<double, std::vector<XData>> sdata;
public:
	template<class TPar, class TL0, class TL1, class DT>
	struct Itr_{
		TPar* parent;
		TL0 level0;
		TL1 level1;

		Itr_& operator++(){
			++level1;
			if (level1 == level0->second.end()){
				++level0;
				if (level0 != parent->sdata.end()){
					level1 = level0->second.begin();
				}
			}
			return *this;
		}
		bool operator!=(const Itr_& it){
			return !(level0 == it.level0 && (level0 == parent->sdata.end() || level1 == it.level1));
		}

		double x() const { return level1->realx; }
		double y() const { return level1->realy; }
		Point point() const { return Point(level1->realx, level1->realy); }
		DT& data(){ return level1->data; }

		//helper static functions
		static Itr_ begin_(TPar* p){
			Itr_ ret;
			ret.parent = p;
			ret.level0 = p->sdata.begin();
			if (p->sdata.size() > 0) ret.level1 = ret.level0->second.begin();
			return ret;
		}
		static Itr_ end_(TPar* p){
			Itr_ ret;
			ret.parent = p;
			ret.level0 = p->sdata.end();
			return ret;
		}
		static Itr_ y_search(TPar* p, TL0 it, double x, double y){
			Itr_ ret;
			ret.parent=p;
			ret.level0 = it;
			assert(it->second.size() > 0);
			for (ret.level1=it->second.begin(); ret.level1!=it->second.end(); ++ret.level1){
				if (fabs(ret.level1->realy-y)<geps && fabs(ret.level1->realx-x)<geps) return ret;
			}
			return end_(p);
		}
		static TL0 find_tl0(TPar* p, double x, TL0 f1){
			if (f1 != p->sdata.end() && f1->first-x < geps) return f1;
			if (f1 != p->sdata.begin()){
				--f1;
				if (x-f1->first<geps) return f1;
			}
			return p->sdata.end();
		}
		static Itr_ find_itr(TPar* p, double x, double y, TL0 f1){
			if (f1 != p->sdata.end() && f1->first-x < 2*geps){
				Itr_ ret = y_search(p, f1, x, y);
				if (ret != end_(p)) return ret;
			}
			if (f1 != p->sdata.begin()){
				--f1;
				if (x-f1->first < 2*geps){
					Itr_ ret = y_search(p, f1, x, y);
					if (ret != end_(p)) return ret;
				}
			}
			return end_(p);
		}
		static Itr_ find_(TPar* p, double x, double y){
			auto f1 = p->sdata.upper_bound(x);
			return find_itr(p, x, y, f1);
		}
		static Itr_ find2_(TPar* p, double x, double y, TL0& fx){
			auto f1 = p->sdata.upper_bound(x);
			fx = find_tl0(p, x, f1);
			return find_itr(p, x, y, f1);
		}
	};

	typedef Itr_<CoordinateMap2D, typename std::map<double, std::vector<XData>>::iterator,
		typename std::vector<XData>::iterator, D> Iterator;
	typedef Itr_<const CoordinateMap2D, typename std::map<double, std::vector<XData>>::const_iterator,
		typename std::vector<XData>::const_iterator, const D> CIterator;
	Iterator begin(){ return Iterator::begin_(this); }
	CIterator begin() const { return CIterator::begin_(this); }
	Iterator end(){ return Iterator::end_(this); }
	CIterator end() const { return CIterator::end_(this); }

	CoordinateMap2D(){ }
	size_t size() const{
		size_t ret=0;
		for (auto it: sdata) ret += it->second.size();
		return ret;
	}

	Iterator find(Point p)  { return find(p.x, p.y); }
	CIterator find(Point p) const { return find(p.x, p.y); }
	Iterator find(double x, double y) { return Iterator::find_(this, x, y);}
	CIterator find(double x, double y) const { return CIterator::find_(this, x, y);}

	//returns if point was added
	bool add(Point key, D val=0) { return add(key.x, key.y, val); }
	bool add(double x, double y, D val=0){
		typename std::map<double, std::vector<XData>>::iterator f1;
		// found equal point
		if (Iterator::find2_(this, x, y, f1) != end()) return false;
		// did not find equal x
		else if (f1 == sdata.end() || fabs(f1->first-x) > geps){
			auto er = sdata.emplace(x, vector<XData>());
			er.first->second.push_back(XData{x, y, val});
			return true;
		} else {
		// found equal x but did not find equal y
			f1->second.push_back(XData{x, y, val});
			return true;
		}
	}

	bool remove(Point key){return remove(key.x, key.y);}
	bool remove(double x, double y){
		Iterator it=find(x, y);
		if (it == end()) return false;
		else{
			remove(it);
			return true;
		}
	}
	void remove(Iterator it){
		it.level0->second.erase(it.level1);
		if (it.level0->second.size() == 0){
			sdata.erase(it.level0);
		}
	}
};

#endif
