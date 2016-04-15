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
	typedef ExtendedSimpleSerialize TSer;
	shared_ptr<TSer> last_ser_result;
	HMCB_SET_PROCNAME("Exporting 3d grid to *.vtk");

	HMCB_SET_DURATION(80, const TSer&, std::string);
	void _run(const TSer& g, std::string fn);

	HMCB_SET_DURATION(
		HMCB_DURATION(ExtendedSimpleSerialize::ConvertExe, Grid) +
		HMCB_DURATION(TGridVTK, TSer, std::string),
		const Grid&, std::string);
	void _run(const Grid& g, std::string fn);
};
extern HMCallback::FunctionWithCallback<TGridVTK> GridVTK;


//signature void BoundaryVTK(const HMGrid3D::Grid& g, std::string* fn);
struct TBoundaryVTK: public HMCallback::ExecutorBase{
	typedef ExtendedSimpleSerialize TSer;
	shared_ptr<TSer> last_ser_result;
	HMCB_SET_PROCNAME("Exporting 3d grid surface to *.vtk");

	HMCB_SET_DURATION(50, TSer, std::string);
	void _run(const TSer&, std::string);


	HMCB_SET_DURATION(
		HMCB_DURATION(HMGrid3D::ExtendedSimpleSerialize::ConvertExe, Grid) +
		HMCB_DURATION(TBoundaryVTK, TSer, std::string),
		const Grid& g, std::string fn
	);
	void _run(const Grid&, std::string);
};
extern HMCallback::FunctionWithCallback<TBoundaryVTK> BoundaryVTK;

//boundary + grid to 2 separate files
//signature void BoundaryVTK(const HMGrid3D::Grid& g, std::string fn);
struct TAllVTK: public HMCallback::ExecutorBase{
	typedef ExtendedSimpleSerialize TSer;
	HMCB_SET_PROCNAME("Exporting 3d data to *.vtk");
	HMCB_SET_DEFAULT_DURATION(
		HMCB_DURATION(ExtendedSimpleSerialize::ConvertExe, Grid)+
		HMCB_DURATION(TBoundaryVTK, TSer, std::string)+ 
		HMCB_DURATION(TGridVTK, TSer, std::string)
	);

	void _run(const Grid& g, std::string fngrid, std::string fnbnd);
};
extern HMCallback::FunctionWithCallback<TAllVTK> AllVTK;

}}

#endif
