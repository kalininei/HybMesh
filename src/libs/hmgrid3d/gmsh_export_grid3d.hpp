#ifndef GMSH_EXPORT_GRID_3D_HPP
#define GMSH_EXPORT_GRID_3D_HPP

#include "hmproject.h"
#include "hmcallback.hpp"
#include "primitives_grid3d.hpp"
#include "serialize_grid3d.hpp"
#include "fluent_export_grid3d.hpp"

namespace HMGrid3D{namespace Export{

//save grid (only tetra, hex, wedge cells
//signature void GridVTK(const HMGrid3D::Grid& g, std::string fn);
struct TGridGMSH: public HMCallback::ExecutorBase{
	HMCB_SET_PROCNAME("Exporting 3d grid to *.gmsh");
	HMCB_SET_DEFAULT_DURATION(100);

	void _run(const SGrid& g, std::string fn);
	void _run(const SGrid& g, std::string fn, BFun func);

};
extern HMCallback::FunctionWithCallback<TGridGMSH> GridGMSH;


}}

#endif


