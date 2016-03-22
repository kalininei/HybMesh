#ifndef HYBMESH_MAPPING_HPP
#define HYBMESH_MAPPING_HPP

#include "grid.h"
#include "hybmesh_contours2d.hpp"
namespace HMGMap{

class MapException: public std::runtime_error{
public:
	MapException(std::string m) noexcept: std::runtime_error(
			std::string("Grid mapping exception: ") + m){};
};

struct Options{
	int fem_nmax;
	int fem_nmin;
	int fem_nrec;
	int fem_nedge;
	std::string snap;   //NO, ADD_VERTICES, SHIFT_VERTICES
	Options(): fem_nmax(100000), fem_nmin(1000), fem_nrec(10000), fem_nedge(3),
	           snap("NO"){}
};

GridGeom MapGrid(const GridGeom& base, const HMCont2D::ECollection& area,
		vector<Point> base_points,
		vector<Point> mapped_points,
		Options opt=Options());


}
#endif
