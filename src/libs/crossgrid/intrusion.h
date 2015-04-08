#ifndef CROSSGRID_INTRUSION_H
#define CROSSGRID_INTRUSION_H

#include "grid.h"
#include "contours.h"


class GridForIntrusion: public GridGeom{
//a hack to access protected GridGeom procedures
//it swaps data (points and cells) with its basic grid
//so the basic grid should not be used before swap_back call
public:
	GridForIntrusion(GridGeom& g){ swap_data(*this, g); }
	void swap_back(GridGeom& g2){ set_indicies(); swap_data(*this, g2);}

	//find a cell containing point
	int find_cell(const Point* p);

	//delete a cell which contains collection of grids
	//adds this collection supporting cells non-overlapping
	//and single connection property
	void intrude(int icell, const vector<GridGeom*>& subs);
	
	//all childs should lie inside one of the parent cells
	static void intrude_grids(GridGeom& parent, const vector<GridGeom*>& childs);
	
	//cover an area between i-th contour and all its child by non-overlapping
	//contours with anti-clockwise direction
	static vector<vector<const Point*>> noncross_filling(const ContoursCollection& col, int i);
};


#endif
