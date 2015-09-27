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
template<class T> using ShpVector = vector<shared_ptr<T>>;


//Macros
#define _DUMMY_FUN_ \
	{\
	printf("DUMMY function:\n"); \
	printf("\t%s\n", __PRETTY_FUNCTION__); \
	printf("\tat %s: %i\n", __FILE__, __LINE__); \
	}

#define _THROW_NOT_IMP_ \
	{\
	printf("function:\t%s\nat\t%s: %i\n", __PRETTY_FUNCTION__, __FILE__, __LINE__); \
	throw std::runtime_error("NOT IMPLEMENTED ERROR");\
	}
#endif
