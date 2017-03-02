#ifndef HMCPORT_C2CPP_HELPER_HPP
#define HMCPORT_C2CPP_HELPER_HPP

#include <vector>
#include <string>
#include <string.h>
#include "bgeom2d.h"
#include <functional>

namespace c2cpp{

template<class A>
void to_ppp(const std::vector<A>& input, int* nret, void*** ret){
	*nret = input.size();
	*ret = new void*[input.size()];
	for (int i=0; i<input.size(); ++i){
		(*ret)[i] = new A(std::move(input[i]));
	}
}

inline void to_char_string(const std::string& str, char** ret){
	*ret = new char[str.size()+10];
	strcpy(*ret, str.c_str());
}

template<class A>
void to_pp(const A& data, void** ret){
	*ret = new A(std::move(data));
}

template<class A>
std::vector<A*> to_pvec(int ndata, void** data){
	std::vector<A*> ret(ndata);
	for (int i=0; i<ndata; ++i){
		ret[i] = static_cast<A*>(data[i]);
	}
	return ret;
}
template<class A>
std::vector<A> to_vec(int ndata, A* data){
	std::vector<A> ret(ndata);
	for (int i=0; i<ndata; ++i){
		ret[i] = data[i];
	}
	return ret;
}

inline vector<Point> to_points2(int np, double* pts){
	vector<Point> ret(np);
	for (int i=0; i<np; ++i){
		ret[i].set(pts[2*i], pts[2*i+1]);
	}
	return ret;
}

inline bool eqstring(const char* _s1, std::string s2){
	if (_s1 == nullptr) return false;
	return s2 == std::string(_s1);
}


inline void assign_boundary_types(int* bnd, int** revdif,
		std::function<void(int, std::map<int, int>&)> whole_assign,
		std::function<int&(int)> bt_by_index){
	std::map<int, int> newbt; //index, newbt
	int sz = *bnd++;
	for (int i=0; i<sz; ++i){
		int nb = *bnd++;
		int btsz = *bnd++;
		for (int j=0; j<btsz; ++j){
			int index = *bnd++;
			if (index < 0) whole_assign(nb, newbt);
			else newbt[index] = nb;
		}
		
	}
	std::map<int, std::vector<int>> rev;
	for (auto kv: newbt){
		int& obt = bt_by_index(kv.first);
		if (obt != kv.second){
			auto fnd = rev.find(kv.first);
			if (fnd == rev.end()){
				fnd = rev.emplace(kv.first, std::vector<int>()).first;
			}
			fnd->second.push_back(obt);
			obt = kv.second;
		}
	}
	//forming return
	int revsz = 1;
	for (auto kv: rev) revsz += 2 + kv.second.size();
	*revdif = new int[revsz];
	int* it = *revdif;
	*it++ = rev.size();
	for (auto kv: rev){
		*it++ = kv.first;
		*it++ = kv.second.size();
		for (int i=0; i<kv.second.size(); ++i){
			*it++ = kv.second[i];
		}
	}
}


}




#endif
