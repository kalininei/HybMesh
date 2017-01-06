#ifndef HYBMESH_TSCALER_HPP
#define HYBMESH_TSCALER_HPP

#include "bgeom2d.h"
#include "bgeom3d.h"
#include "primitives3d.hpp"

namespace Autoscale{

//objects which do scaling to [0, a] square
//and undo it during self destruction
//only objects passed in constructor and add_data() will be automatically unscaled.
class D3{
	HM3D::VertexData _cached;
	ScaleBase3 sc;
	void _do_scale();
	void _undo_scale();
	void _add_to_cache(const HM3D::VertexData& vd);
	double side;
public:
	D3(const D3&) = delete;
	D3(D3&&) = delete;
	D3& operator=(const D3&) = delete;
	D3& operator=(D3&&) = delete;

	//constructors
	D3(const HM3D::VertexData& vd, double a=1.);
	D3(std::initializer_list<const HM3D::VertexData*> vlst, double a=1.);
	D3(const HM3D::FaceData& fd, double a=1.);
	~D3();

	//add to auto unscale list
	void add_data(HM3D::VertexData& vd);
	void add_data(HM3D::FaceData& vd);

	//scale other objects
	void scale(HM3D::VertexData& vd);
	void scale(vector<double>& vd);
	//unscale other objects
	void unscale(HM3D::VertexData& vd);
};


}
#endif
