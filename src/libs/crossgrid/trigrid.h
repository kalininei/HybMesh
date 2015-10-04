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
	ShpVector<Point> ref_points(const vector<double>& dists, double density) const;

	//associates constrain contour with inner contour,
	//places constrain points to contours if necessary
	static vector<HMCont2D::ExtendedTree>
	ConstraintsPreproc(const HMCont2D::ContourTree& cont, 
			const ShpVector<HMCont2D::Contour>& constraints);

	//GModel* gmod
	void FillFromGModel(void* gmod);

public:
	TriGrid(){}
	//constrcut from non-overlapping contours list and additional points
	explicit TriGrid(const ContoursCollection& cont, const vector<double>& lc, double density);

	explicit TriGrid(const HMCont2D::ContourTree& cont, 
			const ShpVector<HMCont2D::Contour>& constraints,
			double h);

	static shared_ptr<TriGrid> FromGmshGeo(const char* fn);

	//procedures get cell centers
	vector<Point> cell_centers() const;

	//laplas smooth
	void smooth(double w);

	static shared_ptr<TriGrid>
	TriangulateArea(const vector<Point>& pts, double h);
	static shared_ptr<TriGrid>
	TriangulateArea(const vector<vector<Point>>& pts, double h);

	static shared_ptr<TriGrid>
	TriangulateAreaConstrained(const vector<vector<Point>>& bnd,
			const vector<vector<Point>>& cns, double h);
};


#endif
