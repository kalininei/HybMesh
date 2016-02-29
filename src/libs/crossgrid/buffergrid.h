#ifndef CROSSGRID_BUFFERGRID_H
#define CROSSGRID_BUFFERGRID_H

#include <map>
#include "grid.h"
#include "hybmesh_contours2d.hpp"

class BufferGrid: public GridGeom{
	//buffer grid which was built within grid 'orig'
	//using buffer 'buffer' from source 'source_cont'
	//it doesn't share points or cells with original grid and source contour
	Contour source_cont;
	double buffer;
	GridGeom* orig;
	vector<const Cell*> orig_cells;  //pointers to cells from orig
	
	void new_edge_points(Edge& e, const vector<double>& wht);

	//filters inner and outer boundary points of the buffer grid
	//std::tuple<
	//        std::set<const Point*>,   //inner_bp - boundary points from source contour
	//        std::set<const Point*>,   //outer_bp - boundary points from internal edges of orig
	//        std::set<const Point*>,   //true_bp  - boundary points from original grid boundary which lies without buffer
	//        std::set<const Point*>    //false_bp - boundary points from original grid boundary which lies whihin buffer
	//> build_bedges() const;

public:
	BufferGrid(GridGeom& main, const PContour& cont, double buffer_size);
	
	int num_orig_cells() const { return orig_cells.size(); }

	//returns contour of the buffer grid with associated points
	//length weights which will be used for triangulation
	std::tuple<
		HMCont2D::ContourTree,
		std::map<Point*, double>
	> boundary_info(bool preserve_true_bp) const;
	
	void update_original() const;
};


#endif
