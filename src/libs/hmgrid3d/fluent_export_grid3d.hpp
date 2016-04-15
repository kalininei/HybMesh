#ifndef FLUENT_EXPORT_GRID3D_HPP
#define FLUENT_EXPORT_GRID3D_HPP
#include "hmproject.h"
#include "primitives_grid3d.hpp"
#include "hmcallback.hpp"

namespace HMGrid3D{namespace Export{
typedef std::function<std::string(int)> BFun;

inline std::string def_bfun(int i){
	return std::string("boundary") + std::to_string(i);
}


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

struct TGridMSH: public HMCallback::ExecutorBase{
	HMCB_SET_PROCNAME("Exporting 3d grid to fluent *.msh");
	HMCB_SET_DEFAULT_DURATION(100);

	void _run(const Grid&, std::string);

	void _run(const Grid&, std::string, BFun);

	HMCB_SET_DURATION(HMCB_DURATION(TGridMSH, Grid, std::string) + 30, 
			Grid, std::string, PeriodicData);
	void _run(const Grid&, std::string, PeriodicData);


	HMCB_SET_DURATION(HMCB_DURATION(TGridMSH, Grid, std::string) + 30, 
			Grid, std::string, BFun, PeriodicData);
	void _run(const Grid&, std::string, BFun, PeriodicData);
};

//instance of TGridMSH for function-like operator() calls
// to use callback call as HMCallback::WithCallback( HMCallback::Fun2, GridMsh, args... );
extern HMCallback::FunctionWithCallback<TGridMSH> GridMSH;



}}
#endif

