#ifndef HMCPORT_C2CPP_HELPER_HPP
#define HMCPORT_C2CPP_HELPER_HPP

#include <vector>
#include <string>
#include <string.h>
#include "bgeom2d.h"

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


}





#endif
