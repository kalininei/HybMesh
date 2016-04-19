#ifndef  REVOLVE_GRID3D_HPP
#define  REVOLVE_GRID3D_HPP

#include "serialize_grid3d.hpp"
#include "grid.h"

namespace HMGrid3D{namespace Constructor{

//returns final grid
HMGrid3D::SGrid RevolveGrid2D(const GridGeom& g2d,
		const vector<double>& phi_coords,
		Point pstart, Point pend, bool is_trian=true,
		std::function<int(int)> side_bt = [](int){ return 1; },
		std::function<int(int)> bt1 = [](int){ return 2; },
		std::function<int(int)> bt2 = [](int){ return 3; });


}}
#endif
