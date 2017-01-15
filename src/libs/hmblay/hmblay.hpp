#ifndef HYBMESH_HMBLAY_HPP
#define HYBMESH_HMBLAY_HPP

#include "primitives2d.hpp"
#include "options.hpp"
#include "hmcallback.hpp"


namespace HMBlay{

class EBuildError: public std::runtime_error{
public:
	EBuildError(std::string m) noexcept: std::runtime_error(
			std::string("Boundary Layer Build Error: ") + m){};
};
HM2D::GridData BuildBLayerGrid(const vector<Input>& opt);

struct TBuildStripeGrid: public HMCallback::ExecutorBase{
	HMCB_SET_PROCNAME("Stripe grid building");
	HMCB_SET_DEFAULT_DURATION(100);

	HM2D::GridData _run(const HM2D::EdgeData& cont,
		const std::vector<double>& partition,
		int tip_algo,
		Point& bl, Point& br, Point& tr, Point& tl);
};
extern HMCallback::FunctionWithCallback<TBuildStripeGrid> BuildStripeGrid;

}

#endif
