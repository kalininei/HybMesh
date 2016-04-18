#ifndef  CONSTRUCT_GRID3D_HPP
#define  CONSTRUCT_GRID3D_HPP

#include "hmgrid3d.hpp"
#include "serialize_grid3d.hpp"
#include "grid.h"

namespace HMGrid3D{namespace Constructor{

//build a cuboid in [0, 0, 0]x[lx, ly, lz] and translate it to leftp.
//boundary types of resulting grid are (for unit cube):
//    x = 0 -> bt = 1
//    x = 1 -> bt = 2
//    y = 0 -> bt = 3
//    y = 1 -> bt = 4
//    z = 0 -> bt = 5
//    z = 1 -> bt = 6
HMGrid3D::Grid Cuboid(HMGrid3D::Vertex leftp, double lx, double ly, double lz, int nx, int ny, int nz);

//sweep xy grid along z vector.
//zcoords represents vector with increasing z coordinate values.
//boundary types:
//   z = zcoords[0] -> bt = 1
//   z = zcoords.back() -> bt = 2
//   all others -> bt= 3
HMGrid3D::Grid SweepGrid2D(const GridGeom& g2d, const vector<double>& zcoords);

//same with supplementary functions defining boundary types
HMGrid3D::Grid SweepGrid2D(const GridGeom& g2d, const vector<double>& zcoords,
		std::function<int(int)> bottom_bt,     //(g2d cell index) -> boundary type
		std::function<int(int)> top_bt,        //(g2d cell index) -> boundary type
		std::function<int(int)> side_bt);      //(g2d edge index) -> boundary type


//returns serialized grid
HMGrid3D::ESS RevolveGrid2D_S(const GridGeom& g2d,
		const vector<double>& phi_coords,
		Point pstart, Point pend, bool is_trian=true,
		std::function<int(int)> side_bt = [](int){ return 1; },
		std::function<int(int)> bt1 = [](int){ return 2; },
		std::function<int(int)> bt2 = [](int){ return 3; });

//returns final grid
HMGrid3D::Grid RevolveGrid2D(const GridGeom& g2d,
		const vector<double>& phi_coords,
		Point pstart, Point pend, bool is_trian=true,
		std::function<int(int)> side_bt = [](int){ return 1; },
		std::function<int(int)> bt1 = [](int){ return 2; },
		std::function<int(int)> bt2 = [](int){ return 3; });

struct Copy{

//deep copies cells, faces, edges.
static HMGrid3D::Grid ShallowVertices(const HMGrid3D::Grid& b);

};



}}

#endif
