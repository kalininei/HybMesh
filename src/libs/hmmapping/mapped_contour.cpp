#include "mapped_contour.hpp"
#include "hmmapping.hpp"

using namespace HMMap::Impl;

double MappedContour::loc2ex_base(double w) const{
	if (ww.size() == 1){
		double ret = w - ww.begin()->first;
		if (ret>1) ret-=1;
		if (ret<0) ret+=1;
		return ret;
	}
	auto it = ww.begin(), itprev = ww.begin();
	double ret = 0;
	while (it!=ww.end() && ret<ww.size()){
		itprev = it++;
		double wprev = itprev->first;
		double wnext = (it == ww.end()) ? ww.begin()->first + 1 : it->first;
		if (it == ww.end() && w < wprev) w+=1;
		if (w>=wprev && w<=wnext){
			ret += (w - wprev)/(wnext - wprev);
			break;
		}
		ret += 1.0;
	}
	assert(ret < ww.size());
	return ret;
}

double MappedContour::loc2ex_mapped(double w) const{
	if (ww_inv.size() == 1){
		double ret = w - ww_inv.begin()->first;
		if (ret>1) ret-=1;
		if (ret<0) ret+=1;
		return ret;
	}
	auto it = ww_inv.begin(), itprev = ww_inv.begin();
	double ret = 0;
	while (it!=ww_inv.end() && ret<ww_inv.size()){
		itprev = it++;
		double wprev = itprev->first;
		double wnext = (it == ww_inv.end()) ? ww_inv.begin()->first + 1 : it->first;
		if (it == ww_inv.end() && w < wprev) w+=1;
		if (w>=wprev && w<=wnext){
			ret += (w - wprev)/(wnext - wprev);
			break;
		}
		ret += 1.0;
	}
	if (ISEQ(ret, ww_inv.size())) ret = 0;
	assert(ret < ww_inv.size());
	return ret;
}

double MappedContour::ex2loc_mapped(double w) const{
	int ibs = std::floor(w);
	auto it = ww.begin();
	for (int i = 0; i<ibs; ++i) ++it;
	double ws = it->second;
	++it;
	if (it == ww.end()) it = ww.begin();
	double we = it->second;
	if (we<=ws) we+=1;
	double t = (w-ibs);
	double ret = (1 - t) * ws + t * we;
	if (ret<0) ret +=1;
	if (ret>1) ret -=1;
	return ret;
}

double MappedContour::ex2loc_base(double w) const{
	int ibs = std::floor(w);
	auto it = ww_inv.begin();
	for (int i = 0; i<ibs; ++i) ++it;
	double ws = it->second;
	++it;
	if (it == ww_inv.end()) it = ww_inv.begin();
	double we = it->second;
	if (we<=ws) we+=1;
	double t = (w-ibs);
	double ret = (1 - t) * ws + t * we;
	if (ret<0) ret +=1;
	if (ret>1) ret -=1;
	return ret;
}

Point MappedContour::map_from_base(Point p) const{
	double w = std::get<1>(base->coord_at(p));
	return map_from_base(w);
}

Point MappedContour::map_from_base(double w) const{
	double ec = loc2ex_base(w);
	double wcont = ex2loc_mapped(ec);
	return HMCont2D::Contour::WeightPoint(*mapped, wcont);
}

Point MappedContour::map_from_mapped(Point p) const{
	double w = std::get<1>(mapped->coord_at(p));
	double ec = loc2ex_mapped(w);
	double wcont = ex2loc_base(ec);
	return HMCont2D::Contour::WeightPoint(*base, wcont);
}

void MappedContour::check_ww() const{
	double w2start = ww.begin()->second;
	double w2prev = -1;
	for (auto wit: ww){
		double w2 = wit.second - w2start;
		if (w2 > 1) w2-=1;
		if (w2 < 0) w2+=1;
		if (w2 <= w2prev) throw HMMap::MapException("Invalid order of points in mapped contour");
		w2prev = w2;
	}
}

void MappedContour::add_connection(Point pbase, Point pmapped){
	double w1 = std::get<1>(base->coord_at(pbase));
	double w2 = std::get<1>(mapped->coord_at(pmapped));
	ww.emplace(w1, w2);
	check_ww();
	ww_inv.emplace(w2, w1);
}

MappedContour* MappedContourCollection::find_by_base(HMCont2D::Contour* bc){
	MappedContour* ret = NULL;
	for (auto s: data){
		if (s->base == bc){ ret = s.get(); break; }
	}
	return ret;
}

MappedContour* MappedContourCollection::find_by_mapped(HMCont2D::Contour* bc){
	MappedContour* ret = NULL;
	for (auto s: data){
		if (s->mapped == bc){ ret = s.get(); break; }
	}
	return ret;
}

MappedContour* MappedContourCollection::insert(HMCont2D::Contour* cbase, HMCont2D::Contour* cmapped){
	auto fnd1 = find_by_base(cbase);
	auto fnd2 = find_by_mapped(cmapped);
	if (fnd1 != 0 && fnd1 == fnd2){
		return fnd1;
	} else if (fnd1 == 0 && fnd2 == 0){
		shared_ptr<MappedContour> n(new MappedContour(cbase, cmapped));
		data.push_back(n);
		return data.back().get();
	} else{
		throw HMMap::MapException("Contour-to-contour links are ambiguous");
	}
}

Point MappedContourCollection::map_from_base(Point p) const{
	//1) find closest contour
	MappedContour* closest=0;
	double closest_w;
	double closest_dist;
	
	for (int i=0; i<data.size(); ++i){
		auto coord = data[i]->base->coord_at(p);
		double dist = std::get<4>(coord);
		if (closest == 0 || dist < closest_dist){
			closest = data[i].get();
			closest_w = std::get<1>(coord);
			closest_dist = dist;
			if (ISZERO(dist)) break;
		}
	}
	
	//2) map
	return closest->map_from_base(closest_w);
}


