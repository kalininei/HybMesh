#ifndef HMGRID3D_MERGE_HPP
#define HMGRID3D_MERGE_HPP

#include "primitives_grid3d.hpp"

namespace HMGrid3D{

//duplicate primitives will be taken from `from` grid
//result will be written to `to` grid.
void MergeGrid(GridData& from, GridData& to,
		const vector<int>& from_vert, const vector<int>& to_vert);


}




#endif
