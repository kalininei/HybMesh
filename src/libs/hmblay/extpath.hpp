#ifndef HYBMESH_HMBLAY_EXTPATH_HPP
#define HYBMESH_HMBLAY_EXTPATH_HPP

#include "options.hpp"

namespace HMBlay{
namespace Impl{


//layer data for point in source path
struct PathPntData{
	PathPntData(Options* o):opt(o){}
	CornerTp tp;
	double angle;
	Vect normal;
	Options* opt;
	Point* p;
	void fill(Point*, Point*, Point*);
	//fill normal, angle according to smoothed contour direction
	//smooth length is calculated acoording to bgrid height
	//if (redefine_types) sets tp = NEGLECTABLE if so
	void set_smooth_angle(int dir, bool redefine_type);
	//same but doesn't make smoothing
	void set_exact_angle(int dir, bool redefine_type);
};

//path with extended info
struct ExtPath: public HM2D::EdgeData{
	typedef HM2D::Edge Ed;

	//PathPnt options for each node. size = ordered_points().size()
	vector<PathPntData> ext_data;

	//contours staring from endpoints and representing boundaries.
	//all angles are between 
	HM2D::EdgeData leftbc, rightbc;

	HM2D::EdgeData* full_source;

	//get methods
	double largest_depth() const;
	int largest_vpart_size() const;


	//makes a partition of path from given lengths of start point to given
	//length of end point. Returns length coordinates of partitions
	vector<double> PathPartition(double, double) const;
	vector<double> VerticalPartition(double) const;

	//assemble from sequence of connected paths given in list of options
	static ExtPath Assemble(const vector<Options*>& src);

	//assembles set of new ExtPath from edges of given one by
	//dividing it by corner type.
	//Corner types remain the same as it were in pth.
	//boundaries are reassembled.
	static vector<ExtPath> DivideByAngle(const ExtPath& pth, CornerTp tp);
	static vector<ExtPath> DivideByAngle(const ExtPath& pth, vector<CornerTp> tps);

	//each resulting section has the same boundary partition
	static vector<ExtPath> DivideByBndPart(const ExtPath& pth);

	//divides at a node which is the closest one to path half point
	static vector<ExtPath> DivideByHalf(const ExtPath& pth);

	//if corner angle section is very short then
	//make it sharp or regular depending on adjecent corner types
	static void ReinterpretCornerTp(ExtPath& pth);

private:
	//builds end perpendicular depending on full source contour
	void FillEndConditions();
	void PerpendicularStart();
	void PerpendicularEnd();

	ExtPath SubPath(const Point* p1, const Point* p2) const;
	ExtPath SubPath(double len1, double len2) const;
	Point* AddPoint(double len);
	HM2D::EdgeData Partition() const;
};









}}


#endif
