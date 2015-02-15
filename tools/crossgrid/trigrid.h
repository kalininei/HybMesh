#ifndef CROSSGRID_TRIGRID_H
#define CROSSGRID_TRIGRID_H

#include <map>
#include "grid.h"
class TriGrid: public GridGeom{
	//properties: edges
	mutable std::set<Edge> _edges;
	std::set<Edge>& edges() const;

	//property: node->node connectivity for inner nodes
	mutable std::map<GridPoint*, vector<GridPoint*>> _nodenodeI;
	std::map<GridPoint*, vector<GridPoint*>>& nodenodeI() const;
	
	//get_refinement points for the meshed area
	shp_vector<Point> ref_points(const vector<double>& dists, double density) const;
public:
	//constrcut from non-overlapping contours list and additional points
	explicit TriGrid(const vector<PContour>& cont, const vector<double>& lc);

	//procedures get cell centers
	vector<Point> cell_centers() const;

	//laplas smooth
	void smooth(double w);
};


#endif
