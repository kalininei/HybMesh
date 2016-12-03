#ifndef HMGRID3D_TETRAHEDRAL_HPP
#define HMGRID3D_TETRAHEDRAL_HPP

#include "hmcallback.hpp"
#include "serialize_grid3d.hpp"
#include "surface_grid3d.hpp"

namespace HMGrid3D{ namespace Mesher{

struct TUnstructuredTetrahedral: public HMCallback::ExecutorBase{
	HMCB_SET_PROCNAME("Tetrahedral meshing");
	HMCB_SET_DEFAULT_DURATION(100);

	GridData _run(const HMGrid3D::Surface& source,
			const HMGrid3D::Surface& sinner,
			const VertexData& pinner, const vector<double>& psizes);

	GridData _run(const HMGrid3D::Surface& source);
	GridData _run(const HMGrid3D::Surface& source,
			const HMGrid3D::Surface& sinner);
	GridData _run(const HMGrid3D::Surface& source,
			const VertexData& pinner, const vector<double>& psizes);
};
extern HMCallback::FunctionWithCallback<TUnstructuredTetrahedral> UnstructuredTetrahedral;



}}

#endif
