#ifndef HYBMESH_HMBLAY_OPTIONS_HPP
#define HYBMESH_HMBLAY_OPTIONS_HPP

#include "primitives2d.hpp"

namespace HMBlay{

enum class BndStepMethod{
	NO_BND_STEPPING,
	CONST_BND_STEP,
	CONST_BND_STEP_KEEP_SHAPE,
	CONST_BND_STEP_KEEP_ALL,
	INCREMENTAL,
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
	HM2D::EdgeData* edges;

	//step along contour length
	double bnd_step;     //for const stepping
	vector<std::pair<Point, double>> bnd_step_basis; //for non-const stepping

	//start, end points of contour tree. Uses hole contour, if they match
	Point start, end;

	//angles
	//0-acute: acute algorithm
	//acute-right: right algorithm
	//right-straight: straight algorithm
	//straight-reentrant: reentrant algorithm
	//reentrant-360: round algorithm
	double acute_angle;
	double right_angle;
	double straight_angle;
	double reentrant_angle;

	//force conformal mapping for all areas.
	//no rectangle approximation if true.
	bool force_conformal;

	//Default values
	Input():
		partition({0}),
		bnd_step_method(BndStepMethod::NO_BND_STEPPING),
		direction(Direction::INNER),
		edges(NULL),
		bnd_step(0.1), bnd_step_basis(),
		start(Point(0,0)),
		end(Point(0, 0)),
		acute_angle(45),
		right_angle(120),
		straight_angle(240),
		reentrant_angle(300),
		force_conformal(false){}

};

namespace Impl{

//corner type: no (start of path), acute, right, straight, reentrant, round,
//neglectable (treat as straight)
enum class CornerTp {NO, ZERO, ACUTE, RIGHT, STRAIGHT, REENTRANT, ROUND, NEGLECTABLE};

//Input options with additional fields. Created by Create From Parent Command.
//All data is deepcopied.
class Options: public Input{
	Options(const Input& other): Input(other){}
	//storage for parent 'edges' member
	shared_ptr<HM2D::EdgeData> __edges_data;

	//these are contours which should be meshed (path)
	//and full contour which contains the path
	//path and full_source are directed in such a way that 
	//grid is build in an inner (left) direction
	//they have deep copied from __edges_data edges and shallow vertices
	HM2D::EdgeData full_source;
	HM2D::EdgeData path;

	//start/end in path
	Point *pnt_start, *pnt_end;

	//all data in [0, 1] square. 'scaling' stores original sizes.
	shared_ptr<ScaleBase> scaling;

	//called to fill pnt_* and path fields
	void Initialize();
public:
	//get
	ScaleBase* get_scaling(){ return scaling.get(); }
	HM2D::EdgeData* get_path(){ return &path; }
	HM2D::EdgeData* get_full_source(){ return &full_source; }

	//treats corner according to this option
	//a is an angle in radians
	CornerTp CornerType(double a);

	//Algorithms
	static vector<Options> CreateFromParent(const vector<Input>& par);
	static vector<vector<Options*>> BuildSequence(vector<Options>& inp);
};

}//Impl




}//HMBlay

#endif
