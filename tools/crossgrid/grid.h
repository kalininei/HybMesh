#ifndef CROSSGRID_GRID_H
#define CROSSGRID_GRID_H

#include "crossgrid.h"
#include "bgeom.h"
#include "contours.h"

class GridPoint: public Point{
	int ind;
public:
	GridPoint(double x, double y, int _ind=0): Point(x,y), ind(_ind){}
	GridPoint(const Point& p): Point(p), ind(0){}
	int get_ind() const { return ind; }

	friend class GridGeom;
};

class Cell{
	vector<GridPoint*> points;
	int ind;
	//reverses points array if neseccary
	void check_ordering();
public:
	explicit Cell(int _ind = 0):ind(_ind){}
	int dim() const { return points.size(); }
	const GridPoint* get_point(int i) const { return points[i]; }
	int get_ind() const { return ind; }
	double area() const;

	friend class GridGeom;
};

struct Edge{
	//nodes indicies
	int p1, p2;
	//cells indicies
	mutable int cell_left, cell_right;
	Edge(int i1=0, int i2=0):
		p1(std::min(i1, i2)), p2(std::max(i1,i2)), 
		cell_left(-1), cell_right(-1){}
	void add_adj_cell(int cell_ind, int i1, int i2) const;
	bool is_boundary() const { return cell_left<0 || cell_right<0; }
};
inline bool operator<(const Edge& e1, const Edge& e2){
	if (e1.p1<e2.p1) return true;
	else if (e1.p1>e2.p1) return false;
	else return (e1.p2<e2.p2);
}

class GridGeom: public Grid{
protected:
	//Data
	shp_vector<GridPoint> points;
	shp_vector<Cell> cells;
	//scaling
	ScaleBase do_scale();
	void do_scale(const ScaleBase& sc);
	void undo_scale(const ScaleBase& sc);
	//contours manipulation
	std::vector<PContour> get_contours() const;
	GridGeom remove_area(const PContour& cont);
	void force_cells_ordering();
	//indexation
	void set_indicies();
	void delete_unused_points();
	//data manipulation
	static void add_point_to_cell(Cell* c, GridPoint* p){ c->points.push_back(p); }
	static void change_point_of_cell(Cell* c, int j, GridPoint* p){ c->points[j] = p; }
	//constructors
	GridGeom(){};
	GridGeom(const GridGeom& g);
	GridGeom& operator=(GridGeom g);
	void add_data(const GridGeom& g);
public:
	//build grid from raw points coordinates array
	//and cells->points connectivity array
	//pts = [x0, y0, x1, y1, ..... ]
	//cls = 
	//	[number of points in cell0, cell0 point0 index, point1 index, ..., 
	//	 number of points in cell1, cell1 point0, cell1 point1,....]
	GridGeom(int Npts, int Ncells, double* pts, int* cls);
	GridGeom(GridGeom&& g);
	~GridGeom(){}

	//number of points
	int n_points() const { return points.size(); }
	//number of cells
	int n_cells() const { return cells.size(); }
	//sum of all cells dimensions
	int n_cellsdim() const;
	
	//data access
	const GridPoint* get_point(int i) const { return points[i].get(); }
	const Cell* get_cell(int i) const { return cells[i].get(); }
	//edges 
	std::set<Edge> get_edges() const;

	//data modify
	//change internal structure of the grid with data from another one.
	//boundary points of this and gg should match
	void change_internal(const GridGeom& gg);

	//static builders
	static GridGeom* cross_grids(GridGeom* gmain, GridGeom* gsec, double buffer_size);
	
	//builds a grid wich is constructed by imposition of gsec onto gmain
	//no bufferzones. 
	//Created grid area equals the intersection of gmain and gsec areas.
	static GridGeom* combine(GridGeom* gmain, GridGeom* gsec);
	static GridGeom* combine2(GridGeom* gmain, GridGeom* gsec);

	friend class BufferGrid;
};














#endif
