#ifndef INCLUDE_HMPROJECT_H
#define INCLUDE_HMPROJECT_H

#include <iostream>
#include <memory>
#include <math.h>
#include <vector>
#include <set>
#include <map>
#include <assert.h>
#include <algorithm>
#include <numeric>
#include <list>

using std::vector;
using std::shared_ptr;
using std::weak_ptr;
template<class T> using ShpVector = vector<shared_ptr<T>>;
template<class T> using WpVector = vector<weak_ptr<T>>;


//Macros
#define _DUMMY_FUN_ \
	{\
	printf("DUMMY function:\n"); \
	printf("\t%s\n", __PRETTY_FUNCTION__); \
	printf("\tat %s: %i\n", __FILE__, __LINE__); \
	}

#define _THROW_NOT_IMP_ \
	{\
	printf("function: %s\nat:       %s: %i\n", __PRETTY_FUNCTION__, __FILE__, __LINE__); \
	throw std::runtime_error("NOT IMPLEMENTED ERROR");\
	}


const double geps = 1e-8;
const double geps2 = 1e-14;
const double gbig = 1.0/geps;
inline bool ISZERO(double v){ return fabs(v)<geps; }
inline bool ISEQ(double a, double b){ return ISZERO(a-b); }
inline bool ISEQLOWER(double x, double y){ return x<y+geps; }
inline bool ISEQGREATER(double x, double y){ return x>y-geps; }
inline bool ISLOWER(double x, double y){ return x<y-geps; }
inline bool ISGREATER(double x, double y){ return x>y+geps; }
inline int SIGN(double x){ if (ISZERO(x)) return 0; else return (x<0) ? -1 : 1; }
inline double sqr(double x){ return x*x; }

//map for gathering coordinates
template<class C>
struct TCoordMap: public std::map<double, C, bool(*)(double, double)>{
	TCoordMap(): std::map<double, C, bool(*)(double, double)>(ISLOWER){}
};
struct TCoordSet: public std::set<double, bool(*)(double, double)>{
	TCoordSet(): std::set<double, bool(*)(double, double)>(ISLOWER){}
};

//positioning constants
const int BOUND = 0;
const int INSIDE = 1;
const int OUTSIDE = -1;

#endif
