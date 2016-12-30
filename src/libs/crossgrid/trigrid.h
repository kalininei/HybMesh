#ifndef CROSSGRID_TRIGRID_H
#define CROSSGRID_TRIGRID_H

#include <map>
#include "grid.h"
#include "hybmesh_contours2d.hpp"

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

	static void CrossesProcessing(
			HMCont2D::ContourTree& cont, 
			ShpVector<HMCont2D::Contour>& constraints,
			std::map<Point*, double>& w,
			HMCont2D::PCollection& apnt,
			double h);

	void guarantee_edges(const vector<HMCont2D::Edge*>& ed);
	void recomb_heal();
public:
	TriGrid(){}
	//constrcut from non-overlapping contours list and additional points
	explicit TriGrid(const ContoursCollection& cont, const vector<double>& lc, double density);

	explicit TriGrid(const HMCont2D::ContourTree& cont, 
			const ShpVector<HMCont2D::Contour>& constraints,
			double h);

	explicit TriGrid(const HMCont2D::ContourTree& cont, 
			const ShpVector<HMCont2D::Contour>& constraints,
			const std::vector<double>& emb_points);

	explicit TriGrid(const HMCont2D::ContourTree& cont, 
			const ShpVector<HMCont2D::Contour>& constraints,
			const std::map<Point*, double>& w, double h);

	static shared_ptr<TriGrid> FromGmshGeo(const char* fn);

	//procedures get cell centers
	vector<Point> cell_centers() const;
	//cell areas
	vector<double> cell_areas() const; 

	//laplas smooth
	void smooth(double w);

	//h is the default size of section
	//w is the function for boundary point -> grid size association
	static shared_ptr<TriGrid>
	TriangulateArea(const vector<Point>& pts, double h);
	static shared_ptr<TriGrid>
	TriangulateArea(const vector<vector<Point>>& pts, double h);
	static shared_ptr<TriGrid>
	TriangulateArea(const HMCont2D::ContourTree& cont, double h);
	static shared_ptr<TriGrid>
	TriangulateArea(const HMCont2D::ContourTree& cont, const std::map<Point*, double>& w, double h);

	static shared_ptr<TriGrid>
	TriangulateAreaConstrained(const vector<vector<Point>>& bnd,
			const vector<vector<Point>>& cns, double h);

	//if h<=0 then h is infinity
	void FillFromTree(
			const HMCont2D::ContourTree& cont, 
			const ShpVector<HMCont2D::Contour>& constraints,
			const vector<Point>& emb_points,
			const std::map<Point*, double>& w,
			double h,
			bool recomb=false);

};

GridGeom QuadGrid(const HMCont2D::ContourTree& cont, 
		const ShpVector<HMCont2D::Contour>& constraints,
		const std::vector<double>& emb_points);
shared_ptr<GridGeom> QuadrangulateArea(const HMCont2D::ContourTree& cont,
		const std::map<Point*, double>& w, double h);


#endif
