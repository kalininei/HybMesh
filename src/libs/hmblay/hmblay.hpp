#ifndef HYBMESH_HMBLAY_HPP
#define HYBMESH_HMBLAY_HPP

#include "grid.h"
#include "options.hpp"


namespace HMBlay{

class EBuildError: public std::runtime_error{
public:
	EBuildError(std::string m) noexcept: std::runtime_error(
			std::string("Boundary Layer Build Error: ") + m){};
};
GridGeom BuildBLayerGrid(const vector<Input>& opt);

}

#endif
