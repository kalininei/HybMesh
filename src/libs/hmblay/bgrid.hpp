#ifndef HYBMESH_BGRID_HPP
#define HYBMESH_BGRID_HPP

#include "grid.h"
#include "hybmesh_contours2d.hpp"
#include "options.hpp"
#include "extpath.hpp"

namespace HMBlay{
namespace Impl{

struct BoundedGrid;

class BGrid: public GridGeom{
protected:
	typedef HMCont2D::Contour Pth;
	typedef HMCont2D::Edge Ed;

	static ExtPath AssembleExtendedPath(vector<Options*>& data);
	static shared_ptr<BGrid> NoSelfIntersections(shared_ptr<BGrid> g);
	static shared_ptr<BGrid> MeshFullPath(const ExtPath& p);
public:
	static shared_ptr<BGrid> MeshSequence(vector<Options*>& data);
	static shared_ptr<BGrid> ImposeBGrids(ShpVector<BGrid>& gg);
};

struct BoundedGrid: public GridGeom{
	virtual vector<GridPoint*> leftside_points() const { _THROW_NOT_IMP_; }
	virtual vector<GridPoint*> rightside_points() const { _THROW_NOT_IMP_; }
};




}//Impl

}//HMBlay



#endif
