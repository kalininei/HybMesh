#ifndef CROSSGRID_BOUNDARY_LAYER_GRID_H
#define CROSSGRID_BOUNDARY_LAYER_GRID_H

#include "grid.h"
#include "hybmesh_contours2d.hpp"


//Represents input data for bl grid building:
//contours, offsets, refinement methods etc
struct BLayerGridInput{

	//perpendicular stepping method
	enum StepMethod {
		STEP0_STEPNUM
	};
	StepMethod step_method;

	//stepping along contours method
	enum BndStepMethod {
		NO_BND_STEPPING,
		CONST_BND_STEP_KEEP_SHAPE,
		CONST_BND_STEP_KEEP_ALL,
		CONST_BND_STEP_NOKEEP,
	};
	BndStepMethod bnd_step_method;

	//input contours tree
	shared_ptr<HMCont2D::ContourTree> tree;
	//direction: INSIDE, OUTSIDE
	int direction;
	//width of grid
	double delta;
	//initial step
	double step0;
	//stepnext = factor*step
	double factor;
	//number of boundary steps
	int step_num;
	//step along contour length
	double bnd_step;
	
	//Constructor with default values
	BLayerGridInput(){
		step_method = STEP0_STEPNUM;
		bnd_step_method = NO_BND_STEPPING;
		direction = INSIDE;
		delta = 1.0;
		step0 = 0.1;
		factor = 1.1;
		step_num = 5;
		bnd_step = 0.1;
	}
};

class BLayerGrid: public GridGeom{
public:
	//throws exception if failed to build a consistent grid
	BLayerGrid(const BLayerGridInput& opt);
};


#endif
