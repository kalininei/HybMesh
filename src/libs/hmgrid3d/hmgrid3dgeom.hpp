#ifndef HMGRID3D_GEOM_HPP
#define HMGRID3D_GEOM_HPP

#include "bgeom2d.h"
#include "primitives_grid3d.hpp"
#include "surface_grid3d.hpp"

namespace HMGrid3D{

struct BoundingBox3D;
struct BoundingSphere3D;

struct BoundingBox3D{
	double xmin, xmax, ymin, ymax, zmin, zmax;

	BoundingBox3D(const Cell& c);
	BoundingBox3D(const Surface& c);
	BoundingBox3D(const VertexData& v);

	//returns INSIDE if bb is fully inside this,
	//        OUTSIDE if fully outside,
	//        BOUND otherwise
	int position(const BoundingBox3D& bb) const;
	int position(const Vertex& v) const;

	//calculate measures.
	//gives 0 if one bb includes the other
	double meas(const BoundingBox3D& bb) const;

	//process
	BoundingBox3D widen(double eps) const;
	BoundingSphere3D inscribe_sphere() const;

private:
	void fill(const ShpVector<Vertex>& v);
};

struct BoundingSphere3D{
	double x, y, z, rad;
	BoundingSphere3D(){}
	BoundingSphere3D(const Surface& c, BoundingBox3D* helper=0);

	//returns INSIDE if bb is fully inside this,
	//        OUTSIDE if fully outside,
	//        BOUND otherwise
	int position(const BoundingBox3D& bb) const;
	int position(const Vertex& v) const;

	//calculate measures.
	//gives 0 if one bb includes the other
	double meas(const BoundingBox3D& bb) const;

	//edit
	BoundingSphere3D widen(double eps) const;

private:
	void fill(const ShpVector<Vertex>& v, BoundingBox3D* helper=0);
	double cmeas(const Vertex& v) const;
};

};


#endif
