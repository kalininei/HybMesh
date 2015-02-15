#ifndef CROSSGRID_BUFFERGRID_H
#define CROSSGRID_BUFFERGRID_H

#include "grid.h"

class BufferGrid: public GridGeom{
	GridGeom* orig;
	vector<const Cell*> orig_cells;
	
	void new_edge_points(Edge& e, const vector<double>& wht);

	//inner and outer boundary points
	//boundary point -> boundary step dictionary
	std::map<const Point*, double> inner_bp;  //boundary points from inner contour
	std::map<const Point*, double> outer_bp;  //boundary points from outer contour
	std::map<const Point*, double> true_bp;   //boundary points from original boundary
	//fill *_bp dictionaries
	void build_bedges(const PContour& cont, double buffer_size);
public:
	BufferGrid(GridGeom& main, const PContour& cont, double buffer_size);
	
	std::tuple<
		std::vector<PContour>,
		std::vector<double>
	> boundary_info() const;
	
	void update_original() const;
};


#endif
