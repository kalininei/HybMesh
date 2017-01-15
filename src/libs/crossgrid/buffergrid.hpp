#ifndef CROSSGRID_BUFFERGRID_H
#define CROSSGRID_BUFFERGRID_H

#include <map>
#include "tree.hpp"

namespace HM2D{ namespace Grid{ namespace Impl{

class BufferGrid: public CellData{
	GridData* orig;
	const EdgeData* source;
	bool preserve_bp;
	double angle0;

	Contour::Tree triangulation_boundary() const;
	void contour_segmentation(EdgeData& cont) const;
	void split_edges(EdgeData&, EdgeData&, EdgeData&, EdgeData&) const;
	VertexData significant_boundary_points(const EdgeData&) const;
	void rebuild_source_edges(EdgeData&, EdgeData&) const;
	EdgeData define_source_edges(const EdgeData&) const;
	int is_on_source(const Point& p) const; //0 - no, 1 - vertex, 2 - on the middle of the edge
	const BoundingBoxFinder& source_finder() const;
	mutable std::unique_ptr<BoundingBoxFinder> _sfinder;
public:
	//extracts cells lying outside const (to the right hand side) not further than buffer_size
	BufferGrid(GridData& main, const EdgeData& source, double buffer_size, bool preserve_bp, double angle0);

	//builds unstructured grid and merges it to original
	//0-triangle, 1-recombined triangle
	void update_original(int filler) const;
};

}}}
#endif
