#ifndef CROSSGRID_BUFFERGRID_H
#define CROSSGRID_BUFFERGRID_H

#include <map>
#include "grid.h"

class BufferGrid: public GridGeom{
	Contour source_cont;
	double buffer;
	GridGeom* orig;
	vector<const Cell*> orig_cells;
	
	void new_edge_points(Edge& e, const vector<double>& wht);

	//inner and outer boundary points
	std::tuple<
		std::set<const Point*>,   //inner_bp - boundary points from inner contour
		std::set<const Point*>,   //outer_bp - boundary points from outer contour
		std::set<const Point*>    //true_bp  - boundary points from original boundary
	> build_bedges() const;

public:
	BufferGrid(GridGeom& main, const PContour& cont, double buffer_size);
	
	int num_orig_cells() const { return orig_cells.size(); }

	std::tuple<
		ContoursCollection,
		std::vector<double>
	> boundary_info(bool preserve_true_bp) const;
	std::tuple<
		ContoursCollection,
		std::vector<double>
	> boundary_info2(bool preserve_true_bp) const;
	
	void update_original() const;
};


#endif
