#ifndef TECPLOT_EXPORT_GRID3D_HPP
#define TECPLOT_EXPORT_GRID3D_HPP
#include "fluent_export_grid3d.hpp"
#include "serialize_grid3d.hpp"
#include "hmcallback.hpp"
namespace HMGrid3D{ namespace Export{

struct TGridTecplot: public HMCallback::ExecutorBase{
	typedef ExtendedSimpleSerialize TSer;
	HMCB_SET_PROCNAME("Exporting 3d grid to tecplot *.dat");
	shared_ptr<TSer> last_ser_result;

	HMCB_SET_DURATION(80, const TSer&, std::string);
	HMCB_SET_DURATION(80, const TSer&, std::string, BFun);
	void _run(const TSer& g, std::string fn, BFun bnd_names=def_bfun);

	HMCB_SET_DURATION(
		HMCB_DURATION(ExtendedSimpleSerialize::ConvertExe, Grid) +
			HMCB_DURATION(TGridTecplot, TSer, std::string),
		const Grid&, std::string);
	HMCB_SET_DURATION(
		HMCB_DURATION(ExtendedSimpleSerialize::ConvertExe, Grid) +
			HMCB_DURATION(TGridTecplot, TSer, std::string),
		const Grid&, std::string, BFun);
	void _run(const Grid& g, std::string fnk, BFun bnd_names=def_bfun);
};

//instance of TGridMSH for function-like operator() calls
// to use callback call as HMCallback::WithCallback( HMCallback::Fun2, GridMsh, args... );
extern HMCallback::FunctionWithCallback<TGridTecplot> GridTecplot;

struct TBoundaryTecplot: public HMCallback::ExecutorBase{
	typedef ExtendedSimpleSerialize TSer;
	HMCB_SET_PROCNAME("Exporting 3d grid surface to tecplot *.dat");
	HMCB_SET_DEFAULT_DURATION(80);

	void _run(const Grid& g, std::string fn, BFun bnd_names=def_bfun);
};

//instance of TGridMSH for function-like operator() calls
// to use callback call as HMCallback::WithCallback( HMCallback::Fun2, GridMsh, args... );
extern HMCallback::FunctionWithCallback<TBoundaryTecplot> BoundaryTecplot;

}}
#endif
