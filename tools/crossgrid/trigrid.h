#ifndef CROSSGRID_TRIGRID_H
#define CROSSGRID_TRIGRID_H

#include <map>
#include "grid.h"
class TriGrid: public GridGeom{
	//initial build procedure
	void build(const vector<PContour>& cont, const shp_vector<Point>& pts);
	
	//properties: edges
	mutable std::set<Edge> _edges;
	std::set<Edge>& edges() const;

	//property: node->node connectivity for inner nodes
	mutable std::map<GPoint*, vector<GPoint*>> _nodenodeI;
	std::map<GPoint*, vector<GPoint*>>& nodenodeI() const;
	
	//get_refinement points for the meshed area
	shp_vector<Point> ref_points(const vector<double>& dists, double density) const;
public:
	//constrcut from non-overlapping contours list and additional points
	explicit TriGrid(const vector<PContour>& cont, const shp_vector<Point>& pts = shp_vector<Point>());
	explicit TriGrid(const vector<PContour>& cont, const vector<Point>& pts);

	//procedures get cell centers
	vector<Point> cell_centers() const;

	//laplas smooth
	void smooth(double w);

	//build better triangle grid within the area of current 
	TriGrid refine_grid(double density);
};


#endif
