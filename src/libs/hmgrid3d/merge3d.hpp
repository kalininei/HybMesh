#ifndef HMGRID3D_MERGE_HPP
#define HMGRID3D_MERGE_HPP

#include "primitives3d.hpp"

namespace HM3D{ namespace Grid{ namespace Algos{

//duplicate primitives will be taken from `from` grid
//result will be written to `to` grid.
void MergeGrid(GridData& from, GridData& to,
		const vector<int>& from_vert, const vector<int>& to_vert);

//deep copied grid, constructed from g1 and g2 will be returned
GridData MergeGrids(const GridData& g1, const GridData& g2);


}}}




#endif
