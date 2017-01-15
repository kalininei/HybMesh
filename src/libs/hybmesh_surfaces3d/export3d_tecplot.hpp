#ifndef TECPLOT_EXPORT_GRID3D_HPP
#define TECPLOT_EXPORT_GRID3D_HPP
#include "export3d_fluent.hpp"
#include "serialize3d.hpp"
#include "hmcallback.hpp"
namespace HM3D{ namespace Export{

struct TGridTecplot: public HMCallback::ExecutorBase{
	HMCB_SET_PROCNAME("Exporting 3d grid to tecplot");
	HMCB_SET_DEFAULT_DURATION(80);

	void _run(const Ser::Grid& g, std::string fn, BFun bnd_names=def_bfun);
	void _run(const GridData& g, std::string fn, BFun bnd_names=def_bfun);

};

//instance of TGridMSH for function-like operator() calls
// to use callback call as HMCallback::WithCallback( HMCallback::Fun2, GridMsh, args... );
extern HMCallback::FunctionWithCallback<TGridTecplot> GridTecplot;

struct TBoundaryTecplot: public HMCallback::ExecutorBase{
	HMCB_SET_PROCNAME("Exporting 3d grid surface to tecplot");
	HMCB_SET_DEFAULT_DURATION(80);

	void _run(const Ser::Grid& g, std::string fn, BFun bnd_names=def_bfun);
};

//instance of TGridMSH for function-like operator() calls
// to use callback call as HMCallback::WithCallback( HMCallback::Fun2, GridMsh, args... );
extern HMCallback::FunctionWithCallback<TBoundaryTecplot> BoundaryTecplot;

}}
#endif
