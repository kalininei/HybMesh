#ifndef HYBMESH_DEBUG_HPP
#define HYBMESH_DEBUG_HPP


#ifndef NDEBUG
#include "hmproject.h"
#include <cxxabi.h>

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

	//print typename of T
	template <class T>
	std::string type_name(){
		typedef typename std::remove_reference<T>::type TR;
		std::unique_ptr<char, void(*)(void*)> own(
			abi::__cxa_demangle(typeid(TR).name(), nullptr, nullptr, nullptr),
			std::free
		);
		std::string r = own != nullptr ? own.get() : typeid(TR).name();
		if (std::is_const<TR>::value) r += " const";
		if (std::is_volatile<TR>::value) r += " volatile";
		if (std::is_lvalue_reference<T>::value) r += "&";
		else if (std::is_rvalue_reference<T>::value) r += "&&";
		return r;
	}
};
#endif
#endif
