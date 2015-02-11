#ifndef CROSSGRID_BUFFERGRID_H
#define CROSSGRID_BUFFERGRID_H

#include "grid.h"

class BufferGrid: public GridGeom{
	GridGeom* orig;
	vector<const Cell*> orig_cells;
public:
	BufferGrid(GridGeom& main, const PContour& cont, double buffer_size);
	
	//modification procedures
	void split_bedges();
	void update_original() const;
};


#endif
