#ifndef HYBMESH_DEBUG_HPP
#define HYBMESH_DEBUG_HPP

#include "hmproject.h"

#ifndef NDEBUG

class HMDebug{
	HMDebug(): ptr(0), i(0) { pre(); }
	void pre();
public:
	static HMDebug& get(){ 
		static HMDebug _dbg;
		return _dbg;
	}

	//printing with tabulation
	static int tabs;
	static void Print(const char* fmt, ...);
	static std::ostream& Cout();

	//ptr, i are free data which can be assign during debugging
	void* ptr;
	int i;
};
#endif
#endif
