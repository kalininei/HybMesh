#ifndef HMGRID3D_TETRAHEDRAL_HPP
#define HMGRID3D_TETRAHEDRAL_HPP

#include "hmcallback.hpp"
#include "surface.hpp"

namespace HM3D{ namespace Mesher{

struct TUnstructuredTetrahedral: public HMCallback::ExecutorBase{
	HMCB_SET_PROCNAME("Tetrahedral meshing");
	HMCB_SET_DEFAULT_DURATION(100);

	GridData _run(const FaceData& source, const FaceData& sinner,
			const VertexData& pinner, const vector<double>& psizes);

	GridData _run(const FaceData& source);
	GridData _run(const FaceData& source, const FaceData& sinner);
	GridData _run(const FaceData& source, const VertexData& pinner,
			const vector<double>& psizes);
};
extern HMCallback::FunctionWithCallback<TUnstructuredTetrahedral> UnstructuredTetrahedral;



}}

#endif
