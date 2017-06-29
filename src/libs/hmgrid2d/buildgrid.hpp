#ifndef HYBMESH_CONSTRUCT_GRIDS_HPP
#define HYBMESH_CONSTRUCT_GRIDS_HPP

#include "primitives2d.hpp"

namespace HM2D{ namespace Grid{ namespace Constructor{


GridData RectGrid01(int Nx, int Ny);
GridData RectGrid(Point p0, Point p1, int Nx, int Ny);
GridData RectGrid(const vector<double>& part_x, const vector<double>& part_y);
//procedures return sequental (x:0->1, y:0->1)
//boundary edges of grids built by RectGrid procedures.
HM2D::EdgeData RectGridBottom(const GridData& gd);
HM2D::EdgeData RectGridRight(const GridData& gd);
HM2D::EdgeData RectGridTop(const GridData& gd);
HM2D::EdgeData RectGridLeft(const GridData& gd);

//rad1 > rad2
GridData Ring(Point p0, double rad1, double rad2, int narc, int nrad);
GridData Circle(Point p0, double rad, int narc, int nrad, bool tri_center);

//uses only those vert which present in vert_cell tabs
GridData FromTab(const VertexData& vert, const vector<vector<int>>& cell_vert);
GridData FromTab(VertexData&& vert, const vector<vector<int>>& cell_vert);

//dim < 2 -> variable cells dimension
//from cell->point table
GridData FromRaw(int npnt, int ncls, double* pnt, int* cls, int dim);

//temporary builds a grid from the collection of cells
//by deleting non-presenting edge-cell connections
//on destroy puts all connections back
class InvokeGrid{
	std::map<int, weak_ptr<Cell>> oldleft;
	std::map<int, weak_ptr<Cell>> oldright;
public:
	InvokeGrid(const InvokeGrid&)=delete;
	InvokeGrid(const CellData& data);
	~InvokeGrid();

	GridData grid;
	void make_permanent();
	
	static GridData Permanent(CellData& data){
		InvokeGrid v(data);
		v.make_permanent();
		return std::move(v.grid);
	}
};

}}}

#endif
