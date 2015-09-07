#ifndef HYBMESH_HMBLAY_OPTIONS_HPP
#define HYBMESH_HMBLAY_OPTIONS_HPP

#include "hybmesh_contours2d.hpp"

namespace HMBlay{

enum class BndStepMethod{
	NO_BND_STEPPING,
	CONST_BND_STEP,
	CONST_BND_STEP_KEEP_SHAPE,
	CONST_BND_STEP_KEEP_ALL,
};
BndStepMethod MethFromString(const char* str);

enum class Direction{
	INNER,
	OUTER,
	BOTH,
};
Direction DirectionFromString(const char* str);

//Represents input data for bl grid building:
//contours, offsets, refinement methods etc
struct Input{
	vector<double> partition;
	//stepping along contours method
	BndStepMethod bnd_step_method;

	//layer direction
	Direction direction;

	//input contours tree. weak reference
	HMCont2D::ECollection* edges;

	//step along contour length
	double bnd_step;

	//rounding off sharp angles
	bool round_off;

	//start, end points of contour tree. Uses hole contour, if they match
	Point start, end;

	//angles
	double sharp_angle;
	double corner_angle;
	double regular_angle;
};

namespace Impl{

//Input options with additional fields. Created by Create From Parent Command.
//All data is deepcopied.
class Options: public Input{
	Options(const Input& other): Input(other){}
	//-- data filled by CreateFromParent and Initialize procedures
	//pool for all edges and points
	shared_ptr<HMCont2D::Container<HMCont2D::ECollection>> __all_data;
	//storage for 'edges' member
	shared_ptr<HMCont2D::ECollection> __edges_data;


	//start, end points and full path for this option.
	//start/end belong to path. Path edges are from __edges_data.
	//equal zero if all contour is used
	Point *pnt_start, *pnt_end;
	HMCont2D::Contour path;

	//all data in [0, 1] square. 'scaling' stores original sizes.
	shared_ptr<ScaleBase> scaling;

	//called to fill pnt_* and path fields
	void Initialize();
public:
	//get
	ScaleBase* get_scaling(){ return scaling.get(); }
	HMCont2D::Contour* get_path(){ return &path; }

	//Algorithms
	static vector<Options> CreateFromParent(const vector<Input>& par);
	static vector<vector<Options*>> BuildSequence(vector<Options>& inp);
};

}//Impl




}//HMBlay

#endif
