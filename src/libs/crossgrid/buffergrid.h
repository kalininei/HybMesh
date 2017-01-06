#ifndef CROSSGRID_BUFFERGRID_H
#define CROSSGRID_BUFFERGRID_H

#include <map>
#include "grid.h"
#include "tree.hpp"

class BufferGrid: public GridGeom{
	//buffer grid which was built within grid 'orig'
	//using buffer 'buffer' from source 'source_cont'
	//it doesn't share points or cells with original grid and source contour
	Contour source_cont;

	double buffer;
	GridGeom* orig;
	vector<const Cell*> orig_cells;  //pointers to cells from orig
	
	void new_edge_points(Edge& e, const vector<double>& wht);
public:
	BufferGrid(GridGeom& main, const PContour& cont, double buffer_size);
	
	int num_orig_cells() const { return orig_cells.size(); }

	//returns contour of the buffer grid with associated points
	//length weights which will be used for triangulation
	std::tuple<
		HM2D::Contour::Tree,
		std::map<Point*, double>
	> boundary_info(bool preserve_true_bp, double angle0) const;
	
	void update_original() const;
};


#endif
