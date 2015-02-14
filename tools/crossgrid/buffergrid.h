#ifndef CROSSGRID_BUFFERGRID_H
#define CROSSGRID_BUFFERGRID_H

#include "grid.h"

class BufferGrid: public GridGeom{
	GridGeom* orig;
	vector<const Cell*> orig_cells;
	
	//buffer grid edges which are boundary for main grid
	struct BndEdge: Edge{
		BndEdge(const Edge& e);
		double w1, w2;
	};
	vector<BndEdge> bedges;

	void build_bedges(const PContour& cont, double buffer_size);
	void new_edge_points(Edge& e, const vector<double>& wht);
public:
	BufferGrid(GridGeom& main, const PContour& cont, double buffer_size);
	
	//modification procedures
	void split_bedges(double density=0.5);

	void update_original() const;
};


#endif
