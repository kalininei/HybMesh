#ifndef EXPORT_GRID_3D_HPP
#define EXPORT_GRID_3D_HPP
#include "hmgrid3d.hpp"
#include "hmcallback.hpp"

namespace HMGrid3D{namespace Export{

class Callback: public HMCallback::Caller2{
	Callback();
public:
	static Callback& instance(std::string s, double duration);
	static Callback& get();

	static void enable(HMCallback::Fun2 f);
	static void disable();
};

template<class Fun, class... Args>
void WithCallback(HMCallback::Fun2 cb, Fun&& f, Args &&... args){
	try{
		Callback::enable(cb);
		f(std::forward<Args>(args)...);
		Callback::disable();
	} catch (...){
		Callback::disable();
		throw;
	}
}

// ===== VTK format
void GridVTK(const HMGrid3D::Grid& g, const char* fn);
void BoundaryVTK(const HMGrid3D::Grid& g, const char* fn);

// ====== Fluent msh format
struct PeriodicDataEntry{
	int bt, bt_shadow;
	Vertex v, v_shadow;
	bool reversed;

	void assemble(HMGrid3D::Grid& g, std::map<Face*, Face*>& outmap);
};
struct PeriodicData{
	//data
	std::vector<PeriodicDataEntry> data;

	//constructing
	PeriodicData(){};

	//* b1 - boundary type for periodic surface
	//* b2 - boundary type for shadow surface
	//* v1/v2 - start points in periodic/shadow surfaces
	//
	//* is_reversed defines connected direction from start point along surface boundary
	//    = true ->
	//      let p1/p2 and p1next/p2next be matched points in periodic/shadow surface contours.
	//      If p1next is located counterclockwise along its contour with respect to p1
	//      then p2next - clockwise with resect to p2, if one looks at corresponding surface so that grid domain lie
	//      beneath it.
	//
	//      is_reversed=true if we want i.g. to match parallel z-surfaces gained by sweep procedure
	//                       so that (x, y, z1) fits (x, y, z2)
	//
	//    = false ->
	//      Otherwise
	void add_condition(int b1, int b2, Vertex v1, Vertex v2, bool is_reversed){
		data.push_back(PeriodicDataEntry{b1, b2, v1, v2, is_reversed});
	}

	//info
	int size() const { return data.size(); }

	//assemble procedure
	HMGrid3D::Grid assemble(const HMGrid3D::Grid& g,
			std::map<Face*, Face*>& outmap);
};

struct GridMSH{
	GridMSH(){};

	//Main functions
	void operator()(const HMGrid3D::Grid& g, const char* fn,
		std::function<std::string(int)> btype_name,
		PeriodicData periodic);

	//Simplified forms:
	//no periodic, all boundaries are called "boundary1", 2, 3 etc.
	void operator()(const HMGrid3D::Grid& g, const char* fn);

	//no periodic
	void operator()(const HMGrid3D::Grid& g, const char* fn,
		std::function<std::string(int)> btype_name);

	//all boundaries are called "boundary1", 2, 3 etc
	void operator()(const HMGrid3D::Grid& g, const char* fn,
		PeriodicData periodic);
};



}}
#endif
