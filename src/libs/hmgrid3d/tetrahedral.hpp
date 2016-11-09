#ifndef HMGRID3D__TETRAHEDRAL_HPP
#define HMGRID3D__TETRAHEDRAL_HPP

#include "serialize_grid3d.hpp"
#include "surface_grid3d.hpp"

namespace HMGrid3D{ namespace Mesher {

SGrid UnstructedTetrahedral(const ShpVector<HMGrid3D::Surface>& srfs);



}}

#endif
