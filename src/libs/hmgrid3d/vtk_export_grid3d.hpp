#ifndef VTK_EXPORT_GRID_3D_HPP
#define VTK_EXPORT_GRID_3D_HPP

#include "hmproject.h"
#include "hmcallback.hpp"
#include "primitives_grid3d.hpp"
#include "serialize_grid3d.hpp"

namespace HMGrid3D{namespace Export{

//save grid (only tetra, hex, wedge cells
//signature void GridVTK(const HMGrid3D::Grid& g, std::string fn);
struct TGridVTK: public HMCallback::ExecutorBase{
	HMCB_SET_PROCNAME("Exporting 3d grid to *.vtk");

	HMCB_SET_DURATION(80, const SGrid&, std::string);
	void _run(const SGrid& g, std::string fn);

};
extern HMCallback::FunctionWithCallback<TGridVTK> GridVTK;


//signature void BoundaryVTK(const HMGrid3D::Grid& g, std::string* fn);
struct TBoundaryVTK: public HMCallback::ExecutorBase{
	HMCB_SET_PROCNAME("Exporting 3d grid surface to *.vtk");

	HMCB_SET_DURATION(50, SGrid, std::string);
	void _run(const SGrid&, std::string);

};
extern HMCallback::FunctionWithCallback<TBoundaryVTK> BoundaryVTK;

//boundary + grid to 2 separate files
struct TAllVTK: public HMCallback::ExecutorBase{
	HMCB_SET_PROCNAME("Exporting 3d data to *.vtk");
	HMCB_SET_DEFAULT_DURATION(
		HMCB_DURATION(TBoundaryVTK, SGrid, std::string)+ 
		HMCB_DURATION(TGridVTK, SGrid, std::string)
	);

	void _run(const SGrid& g, std::string fngrid, std::string fnbnd);
};
extern HMCallback::FunctionWithCallback<TAllVTK> AllVTK;

}}

#endif
