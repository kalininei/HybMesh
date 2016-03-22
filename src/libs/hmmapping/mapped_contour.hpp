#ifndef MAPPED_CONTOUR_HPP
#define MAPPED_CONTOUR_HPP
#include "hybmesh_contours2d.hpp"

namespace HMGMap{
namespace Impl{

class MappedContourCollection;

class MappedContour{
	const HMCont2D::Contour* base;
	const HMCont2D::Contour* mapped;
	std::map<double, double> ww;
	void check_ww() const;
	//extended coordinate = int(index of base point) + float(edge coordinate between base points)
	double loc2ex_base(double w) const;
	double ex2loc_mapped(double w) const;
public:
	MappedContour(const HMCont2D::Contour* c1, const HMCont2D::Contour* c2): base(c1), mapped(c2){};

	Point map_from_base(Point p) const;
	void add_connection(Point pbase, Point pmapped);
	friend class MappedContourCollection;
};


class MappedContourCollection: public ShpVector<MappedContour>{
	ShpVector<MappedContour> data;
public:
	int entry_num() const { return data.size(); }
	MappedContour* find_by_base(HMCont2D::Contour* bc);
	//inserts iff cbase->cmapped pair was not found.
	MappedContour* insert(HMCont2D::Contour* cbase, HMCont2D::Contour* cmapped);
};





}}
#endif
