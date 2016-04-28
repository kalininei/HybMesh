#ifndef MAPPED_CONTOUR_HPP
#define MAPPED_CONTOUR_HPP
#include "hybmesh_contours2d.hpp"

namespace HMGMap{
namespace Impl{

class MappedContourCollection;

class MappedContour{
	const HMCont2D::Contour* base;
	const HMCont2D::Contour* mapped;
	std::map<double, double> ww;  //base->mapped
	std::map<double, double> ww_inv; //mapped->base
	void check_ww() const;
	//extended coordinate = int(index of base point) + float(edge coordinate between base points)
	double loc2ex_base(double w) const;
	double loc2ex_mapped(double w) const;
	double ex2loc_base(double w) const;
	double ex2loc_mapped(double w) const;
public:
	MappedContour(const HMCont2D::Contour* c1, const HMCont2D::Contour* c2): base(c1), mapped(c2){};

	const HMCont2D::Contour* get_base() const { return base; }
	const HMCont2D::Contour* get_mapped() const { return mapped; }

	Point map_from_base(Point p) const; //returns point on mapped contour
	Point map_from_mapped(Point p) const;  //returns point on base contour
	void add_connection(Point pbase, Point pmapped);
	friend class MappedContourCollection;
};


class MappedContourCollection{
	ShpVector<MappedContour> data;
public:
	int entry_num() const { return data.size(); }
	MappedContour* find_by_base(HMCont2D::Contour* bc);
	MappedContour* find_by_mapped(HMCont2D::Contour* bc);
	const MappedContour* get(int i) const { return data[i].get(); }
	//inserts iff cbase->cmapped pair was not found.
	MappedContour* insert(HMCont2D::Contour* cbase, HMCont2D::Contour* cmapped);
};





}}
#endif
