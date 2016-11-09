#ifndef HMGRID3D_UNITE_GRID3D_HPP
#define HMGRID3D_UNITE_GRID3D_HPP
#include "serialize_grid3d.hpp"

namespace HMGrid3D{ namespace Unite{

SGrid Superimpose(const SGrid& base, const SGrid& imp, double buf);

}}


#endif
