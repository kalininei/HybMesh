#ifndef FLUENT_EXPORT_GRID3D_HPP
#define FLUENT_EXPORT_GRID3D_HPP
#include "hmproject.h"
#include "hmcallback.hpp"
#include "serialize3d.hpp"
#include "treverter3d.hpp"

namespace HM3D{namespace Export{

typedef std::function<std::string(int)> BFun;

inline std::string def_bfun(int i){
	return std::string("boundary") + std::to_string(i);
}


// ====== Fluent msh format
struct PeriodicDataEntry{
	int bt, bt_shadow;
	Vertex v, v_shadow;
	bool reversed;

	void assemble(GridData& g, std::map<Face*, Face*>& outmap);
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
	GridData assemble(const HM3D::GridData& g,
			std::map<Face*, Face*>& outmap);
};

struct TGridMSH: public HMCallback::ExecutorBase{
	HMCB_SET_PROCNAME("Exporting 3d grid to fluent file");
	HMCB_SET_DEFAULT_DURATION(100);

	void _run(const Ser::Grid&, std::string);
	void _run(const Ser::Grid&, std::string, BFun);

	HMCB_SET_DURATION(HMCB_DURATION(TGridMSH, Ser::Grid, std::string) + 30, 
			Ser::Grid, std::string, PeriodicData);
	void _run(const Ser::Grid&, std::string, PeriodicData);


	HMCB_SET_DURATION(HMCB_DURATION(TGridMSH, Ser::Grid, std::string) + 30, 
			const Ser::Grid&, std::string, BFun, PeriodicData);
	void _run(const Ser::Grid&, std::string, BFun, PeriodicData);

	// ===== GridData versions
	void _run(const GridData& a, std::string b){ return _run(Ser::Grid(a), b); }
	void _run(const GridData& a, std::string b, BFun c) { return _run(Ser::Grid(a), b, c); }

	HMCB_SET_DURATION(HMCB_DURATION(TGridMSH, Ser::Grid, std::string) + 30, 
			GridData, std::string, PeriodicData);
	void _run(const GridData& a, std::string b, PeriodicData c){ return _run(Ser::Grid(a), b, c); }


	HMCB_SET_DURATION(HMCB_DURATION(TGridMSH, Ser::Grid, std::string) + 30, 
			GridData, std::string, BFun, PeriodicData);
	void _run(const GridData& a, std::string b, BFun c, PeriodicData d){ return _run(Ser::Grid(a), b, c, d); }

};

//instance of TGridMSH for function-like operator() calls
// to use callback call as HMCallback::WithCallback( HMCallback::Fun2, GridMsh, args... );
extern HMCallback::FunctionWithCallback<TGridMSH> GridMSH;



}}
#endif

