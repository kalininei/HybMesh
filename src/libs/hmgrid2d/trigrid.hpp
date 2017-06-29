#ifndef HYBMESH_TRIGRID_HPP
#define HYBMESH_TRIGRID_HPP
#include "hmcallback.hpp"
#include "primitives2d.hpp"
#include "contour_tree.hpp"
#include "nodes_compare.h"

namespace HM2D{ namespace Mesher{

//makes smooth 1d segmentation of source keeping all points.
//creates new edges, but uses source vertices if possible
Contour::Tree PrepareSource(const Contour::Tree& source, double defsize = -1);
Contour::Tree PrepareSource(const EdgeData& source, double defsize = -1);

//does calculations for a single closed contour with given source points
//!!!! all edges and vertices with id = 1 will not be changed and used as size sources.
//creates new edges, but uses source vertices if possible
EdgeData RepartSourceById(const EdgeData& source, const vector<std::pair<Point, double>>& src);

//unstructured meshing procedures fills domain using existing boundary segmentation.
//all detached tree nodes will be treated as constraints
struct TUnstructuredTriangle: public HMCallback::ExecutorBase{
	HMCB_SET_PROCNAME("Triangulation");
	HMCB_SET_DEFAULT_DURATION(100);

	GridData _run(const Contour::Tree& source);
	GridData _run(const Contour::Tree& source, const CoordinateMap2D<double>& embedded);
};
extern HMCallback::FunctionWithCallback<TUnstructuredTriangle> UnstructuredTriangle;

struct TUnstructuredTriangleRecomb: public HMCallback::ExecutorBase{
	HMCB_SET_PROCNAME("Triangulation with recombine");
	HMCB_SET_DEFAULT_DURATION(110);

	GridData _run(const Contour::Tree& source);
	GridData _run(const Contour::Tree& source, const CoordinateMap2D<double>& embedded);
};
extern HMCallback::FunctionWithCallback<TUnstructuredTriangleRecomb> UnstructuredTriangleRecomb;


}}
#endif
