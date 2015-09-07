#ifndef HYBMESH_BGRID_HPP
#define HYBMESH_BGRID_HPP

#include "grid.h"
#include "hybmesh_contours2d.hpp"
#include "options.hpp"

namespace HMBlay{
namespace Impl{


//corner type: no (start of path), sharp, corner (right), regular (plain), obtuse
enum class CornerTp {NO, SHARP, CORNER, REGULAR, OBTUSE};

//layer data for point in source path
struct PathPntData{
	CornerTp tp;
	Vect normal;
	Options* opt;
	double angle; //angle between previous/next edge
};
PathPntData BuildPathPntData(Point* p1, Point* p2, Point* p3, Options* opt);

//path with extended info
struct ExtPath: public HMCont2D::Contour{
	typedef HMCont2D::Edge Ed;

	//PathPnt options for each node. size = ordered_points().size()
	vector<PathPntData> ext_data;

	//assemble from sequence of connected paths given in list of options
	static ExtPath Assemble(const vector<Options*>& src);
};

class BGrid: public GridGeom{
protected:
	typedef HMCont2D::Contour Pth;
	typedef HMCont2D::Edge Ed;
	typedef std::map<ExtPath*, vector<ExtPath>> TEpMap;
	//--- data
	//cell priority:
	//cells which are close to source have high priority.
	//the closest have priority = HIGHEST_PRIORIRY
	const int HIGHEST_PRIORITY = 1000;
	std::map<Cell*, int> priority;

	//--- assembling pathes
	static ExtPath AssembleExtendedPath(vector<Options*>& data);
	//--- dividing pathes
	static TEpMap RemoveObtuse(vector<ExtPath*>&&);
	static TEpMap RemoveSharp(vector<ExtPath*>&&);
	static TEpMap RemoveCorner(vector<ExtPath*>&&);
	//--- merging grids
	static BGrid JoinCornerNodes(vector<std::pair<ExtPath*, BGrid*>>& inp);
	static BGrid JoinSharpNodes(vector<std::pair<ExtPath*, BGrid*>>& inp);
	static BGrid JoinObtuseNodes(vector<std::pair<ExtPath*, BGrid*>>& inp);
public:
	static ShpVector<BGrid> MeshSequence(vector<Options*>& data);
	static BGrid ImposeBGrids(ShpVector<BGrid>& gg);
};


//grid with rectangular structure
class SimpleBGrid: public BGrid{
	//each stencil entry is a line, perpendicular to source
	vector<vector<GridPoint*>> stencil;
	bool is_closed;
public:
	SimpleBGrid(ExtPath& pth);
};


}//Impl

}//HMBlay



#endif
