#ifndef CROSSGRID_BOUNDARY_LAYER_GRID_H
#define CROSSGRID_BOUNDARY_LAYER_GRID_H

#include "grid.h"
#include "hybmesh_contours2d.hpp"


//Represents input data for bl grid building:
//contours, offsets, refinement methods etc
struct BLayerGridInput{

	vector<double> partition;
	//stepping along contours method
	enum BndStepMethod {
		NO_BND_STEPPING,
		CONST_BND_STEP,
		CONST_BND_STEP_KEEP_SHAPE,
		CONST_BND_STEP_KEEP_ALL,
	};
	BndStepMethod bnd_step_method;

	//input contours tree. weak reference
	HMCont2D::ContourTree* tree;
	//direction: INSIDE, OUTSIDE
	int direction;

	//step along contour length
	double bnd_step;
	//rounding off sharp angles
	bool round_off;
	//minimal sharp angle (radians)
	double minsharp;

	//Constructor with default values
	BLayerGridInput(){
		bnd_step_method = NO_BND_STEPPING;
		direction = INSIDE;
		bnd_step = 0.1;
	}
};

class BLayerGrid: public GridGeom{
public:
	//throws exception if failed to build a consistent grid
	BLayerGrid(const BLayerGridInput& opt);
};


#endif
