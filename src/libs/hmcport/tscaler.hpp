#ifndef HYBMESH_TSCALER_HPP
#define HYBMESH_TSCALER_HPP

#include "bgeom2d.h"
#include "bgeom3d.h"
#include "primitives3d.hpp"
#include "primitives2d.hpp"

namespace Autoscale{

//objects which do scaling to [0, a] square
//and undo it during self destruction
//only objects passed in constructor and add_data() will be automatically unscaled.
class D2{
	vector<HM2D::EdgeData*> _einp;
	vector<HM2D::GridData*> _ginp;
	ScaleBase sc;
public:
	//no defaults just in case to prevent multiple unscaling
	//on destruction.
	D2(const D2&) = delete;
	D2(D2&&) = delete;
	D2& operator=(const D2&) = delete;
	D2& operator=(D2&&) = delete;

	~D2();
	//simple constructors
	D2(vector<Point>& p, double a=1.);
	//constructors with automatic unscaling on destruction
	D2(HM2D::EdgeData* vd, double a=1.);
	D2(HM2D::GridData* vd, double a=1.);
	D2(const vector<HM2D::EdgeData*>& vd, double a=1.);

	//scale with auto unscale
	void add_data(HM2D::EdgeData* vd);
	void add_data(HM2D::GridData* vd);
	void add_data(const vector<HM2D::EdgeData*>& vd);

	//simple scale other objects
	void scale(vector<Point>& pts);
	void scale(vector<double>& lens);
	void scale(Point& p);
	void scale(double& a);
	
	//unscale other objects
	void unscale(HM2D::EdgeData* vd);
	void unscale(HM2D::GridData* vd);
	void unscale(HM3D::GridData* vd);

	//misc
	ScaleBase get_scale() const{ return sc; }
};

class D3{
	vector<HM3D::GridData*> _ginp;
	vector<HM3D::FaceData*> _finp;
	ScaleBase3 sc;
public:
	D3(const D3&) = delete;
	D3(D3&&) = delete;
	D3& operator=(const D3&) = delete;
	D3& operator=(D3&&) = delete;

	~D3();
	//constructors with automatic unscaling on destruction
	D3(HM3D::FaceData* vd, double a=1.);
	D3(HM3D::GridData* vd, double a=1.);
	D3(const vector<HM3D::GridData*>& vd, double a=1.);
	D3(const vector<HM3D::FaceData*>& vd, double a=1.);

	//scale other objects
	//...
	//unscale other objects
	void unscale(HM3D::GridData* g);

};

}
#endif
