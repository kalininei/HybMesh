#ifndef CROSSGRID_PROCGRID_H
#define CROSSGRID_PROCGRID_H
//This is a new style interface.
//Old interface should be completely changed to new one during grid refactoring.
#include "hybmesh_contours2d.hpp"
#include <sstream>
#include "grid.h"
namespace GGeom{

// === exceptions
//basic exception
class EException: public std::runtime_error{
public:
	EException(std::string m) noexcept: std::runtime_error(
			std::string("Grid geometry exception: ") + m){};
};
//point is out of area
class EOutOfArea: public EException{
	static std::string Msg(Point p){
		std::ostringstream s;
		s<<"Point "<<p<<" lies outside grid area"<<std::endl;
		return s.str();
	};
public:
	EOutOfArea(Point p) noexcept: EException(Msg(p)){};
};


// === builders
struct Constructor{

static GridGeom EmptyGrid();
static GridGeom RectGrid(const vector<double>& part_x, vector<double>& part_y);
static GridGeom RectGrid(Point p0, Point p1, int Nx, int Ny);
//Nx, Ny - number of nodes in x, y directions
static GridGeom RectGrid01(int Nx, int Ny);
//rad1 > rad2
static GridGeom Ring(Point p0, double rad1, double rad2, int narc, int nrad);
static GridGeom Circle(Point p0, double rad, int narc, int nrad, bool tri_center);
static GridGeom DeepCopy(const GridGeom& g);   //g should be indexed
//builds grid from given cells indicies.
//policy:
// 0 - make a deep copy of points and cells
// 1 - shallow copy + set indiceis from returned grid
// 2 - shallow copy + leave indicies untouched
//works even if g is not well indexed
static GridGeom ExtractCells(const GridGeom& g, const std::vector<int>& cind, int policy=0);
};


// === modifiers
struct Modify{

static void RemoveCells(GridGeom& grid, const std::vector<const Cell*>& cls);
static void AddCell(GridGeom& grid, const std::vector<Point>& cell);
//primitives modifications
static void PointModify(GridGeom& grid, std::function<void(GridPoint*)> fun);
static void CellModify(GridGeom& grid, std::function<void(Cell*)> fun);
static void ClearAll(GridGeom& grid);
//adds data
//!! as a result of ShallowAdd points and cells indicies are renumbered according to
//'to' grid index. Be sure that 'from' grid is not used after this procedure.
static void ShallowAdd(const GridGeom* from, GridGeom* to);
static void ShallowAdd(const ShpVector<GridPoint>& from, GridGeom* to);
//adds certain cells indicies. if ind_to then sets indicies from 'to' grid. Else doesn't change them
static void ShallowAdd(const GridGeom* from, GridGeom* to, const std::vector<int>& icells,
		bool ind_to=true);
static void DeepAdd(const GridGeom* from, GridGeom* to);

//Set of single connected meshes. Cells and Points are shallow copies
//of the original. All indicies left unchanged
static vector<GridGeom> SubGrids(const GridGeom& grid);

//if boundary edge points (p1, p2) lie on cont then
//all significant contour points between (p1, p2) will present in grid
//grid should be located to left of contour
//snap_nodes is a list of points which will be snapped to contour before procedure starts
static void SnapToContour(GridGeom& grid, const HMCont2D::Contour& cont,
		const std::vector<GridPoint*>& snap_nodes);
//shifts boundary grid node to significant contour vertex
//if it is non-significant by itself. Otherwise does nothing
static void ShiftToContour(GridGeom& grid, const HMCont2D::Contour& cont,
		const std::vector<GridPoint*>& snap_nodes);

//no complicated boundary cell edges
//angle is between [0, M_PI]
//0 -- delete only non-significant edge points
//pi -- delete all intermediate points
static void SimplifyBoundary(GridGeom& grid, double angle);

private:
	struct _ShiftSnapPreCalc;
};

// === grid structure information
struct Info{

//points
static ShpVector<GridPoint> SharePoints(const GridGeom& grid);
static ShpVector<GridPoint> SharePoints(const GridGeom& grid, const vector<int>& indicies);
static ShpVector<GridPoint> BoundaryPoints(const GridGeom& grid);
//collection of all outer contours
static HMCont2D::ContourTree Contour(const GridGeom& grid);
//only the first contour of tree.
//Used for grids which are definitly singly connected
static HMCont2D::Contour Contour1(const GridGeom& grid);
//contour from grid cell
static HMCont2D::Contour CellContour(const GridGeom& grid, int cell_ind);
//Build a bounding box
static BoundingBox BBox(const GridGeom& grid, double eps=geps);
//calculate skewness
static vector<double> Skewness(const GridGeom& grid);
//calculates area as the sum of all cells areas
static double Area(const GridGeom& grid);
//gets cells areas as vector
static vector<double> CellAreas(const GridGeom& grid);
//checks all cells for correct direction and non-intersection
static bool Check(const GridGeom& grid);

//Finders
class CellFinder{
	const GridGeom* grid;
	const int Nx, Ny;
	double Hx, Hy;
	Point p0;
	vector<vector<const Cell*>> cells_by_square;

	int GetSquare(const Point& p) const;
	std::set<int> IndSet(const BoundingBox& bbox) const;
public:
	CellFinder(const GridGeom* g, int nx, int ny);
	//Find cell containing point.
	//Throws EOutOfArea if failed to find the point.
	const Cell* Find(const Point& p) const;
	const Cell* FindExcept(const Point& p, const std::set<const Cell*>& exc) const;

	//vector of cells which can contain given point
	const vector<const Cell*>& CellCandidates(const Point& p) const;
	const vector<const Cell*>& CellsBySquare(int i) const { return cells_by_square[i]; }
	//adds cell to info vector
	void AddCell(const Cell* c);
};

};

struct Repair{

//merges congruent, deletes unused, forces cells rotation
//edge directing, merging points, no unused points
static void Heal(GridGeom& grid);
static bool HasSelfIntersections(const GridGeom& grid);

};




}

#endif
