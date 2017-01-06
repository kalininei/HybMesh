#ifndef MAPPED_CONTOUR_HPP
#define MAPPED_CONTOUR_HPP
#include "primitives2d.hpp"

namespace HMMap{
namespace Impl{

class MappedContourCollection;

class MappedContour{
	const HM2D::EdgeData* base;
	const HM2D::EdgeData* mapped;
	std::map<double, double> ww;  //base->mapped
	std::map<double, double> ww_inv; //mapped->base
	void check_ww() const;
	//extended coordinate = int(index of base point) + float(edge coordinate between base points)
	double loc2ex_base(double w) const;
	double loc2ex_mapped(double w) const;
	double ex2loc_base(double w) const;
	double ex2loc_mapped(double w) const;

	HM2D::EdgeData rev_mapped;
	const HM2D::EdgeData* orig_mapped;
	Point map_from_base(double w) const; //returns point on mapped contour from base weight 
public:
	MappedContour(const HM2D::EdgeData* c1, const HM2D::EdgeData* c2, bool reversed);

	const HM2D::EdgeData* get_base() const { return base; }
	const HM2D::EdgeData* get_mapped() const { return mapped; }
	const HM2D::EdgeData* get_orig_mapped() const { return orig_mapped; }

	Point map_from_base(Point p) const; //returns point on mapped contour

	Point map_from_mapped(Point p) const;  //returns point on base contour
	void add_connection(Point pbase, Point pmapped);
	friend class MappedContourCollection;
};


class MappedContourCollection{
	ShpVector<MappedContour> data;
	const bool reversed;
public:
	MappedContourCollection(bool reversed): reversed(reversed){}
	int entry_num() const { return data.size(); }
	bool is_reversed() const {return reversed; }
	MappedContour* find_by_base(HM2D::EdgeData* bc);
	MappedContour* find_by_mapped(HM2D::EdgeData* bc);
	const MappedContour* get(int i) const { return data[i].get(); }
	//inserts iff cbase->cmapped pair was not found.
	MappedContour* insert(HM2D::EdgeData* cbase, HM2D::EdgeData* cmapped);

	Point map_from_base(Point p) const;
};





}}
#endif
