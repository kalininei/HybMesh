#ifndef GMSH_EXPORT_GRID_3D_HPP
#define GMSH_EXPORT_GRID_3D_HPP

#include "hmproject.h"
#include "hmcallback.hpp"
#include "primitives3d.hpp"
#include "serialize3d.hpp"
#include "export3d_fluent.hpp"

namespace HM3D{namespace Export{

//save grid (only tetra, hex, wedge cells
//signature void GridVTK(const HMGrid3D::Grid& g, std::string fn);
struct TGridGMSH: public HMCallback::ExecutorBase{
	HMCB_SET_PROCNAME("Exporting 3d grid to gmsh");
	HMCB_SET_DEFAULT_DURATION(100);

	void _run(const Ser::Grid& g, std::string fn);
	void _run(const Ser::Grid& g, std::string fn, BFun func);

};
extern HMCallback::FunctionWithCallback<TGridGMSH> GridGMSH;


}}

#endif


