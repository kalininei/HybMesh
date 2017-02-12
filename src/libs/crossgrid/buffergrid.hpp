#ifndef CROSSGRID_BUFFERGRID_H
#define CROSSGRID_BUFFERGRID_H

#include <map>
#include "contour_tree.hpp"

namespace HM2D{ namespace Grid{ namespace Impl{

class BufferGrid: public CellData{
	GridData* orig;
	const EdgeData* source;
	bool preserve_bp;
	double angle0;

	//orig points lying on source which were deleted
	VertexData badpoints;

	vector<EdgeData> build_contacts(const EdgeData& bedges) const;
	Contour::Tree offset_contact(const Contour::Tree& outer_bnd, const EdgeData& contact, double bsize) const;
	Contour::Tree offset_closed_contact(const Contour::Tree& outer_bnd, const EdgeData& contact, double bsize) const;
	Contour::Tree offset_open_contact(const Contour::Tree& outer_bnd, const EdgeData& contact, double bsize) const;
	void remove_vertices_from_orig(const VertexData& badp);
	Contour::Tree build_buffer_zone(const Contour::Tree& outer_bnd, double buffer_size) const;
	Contour::Tree triangulation_boundary();

	struct splitR{
		EdgeData bsource, bgrid, bboundary;
		VertexData keep_points;
	};
	splitR split_edges(EdgeData&);
	void build_size_sources(const splitR&, vector<std::pair<Point, double>>&);
	void contour_segmentation(EdgeData& cont, const splitR& splitted,
			const vector<std::pair<Point, double>>& src);
	VertexData significant_boundary_points(const EdgeData&) const;
	void rebuild_source_edges(EdgeData&, EdgeData&);
	EdgeData define_source_edges(const EdgeData&) const;
	int is_on_source(const Point& p) const; //0 - no, 1 - vertex, 2 - somewhere on the edge
	const BoundingBoxFinder& source_finder() const;
	mutable std::unique_ptr<BoundingBoxFinder> _sfinder;
public:
	//extracts cells lying outside const (to the right hand side) not further than buffer_size
	BufferGrid(GridData& main, const EdgeData& source, double buffer_size, bool preserve_bp, double angle0);

	//builds unstructured grid and merges it to original
	//0-triangle, 1-recombined triangle
	void update_original(int filler);
};

}}}
#endif
