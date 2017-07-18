#ifndef HYBMESH_INSCRIBE_GRID_HPP
#define HYBMESH_INSCRIBE_GRID_HPP

#include "primitives2d.hpp"
#include "hmcallback.hpp"
#include "contour_tree.hpp"

namespace HM2D{ namespace Grid{ namespace Algos{

enum class SubstractCellsAlgo{
	               //Which cell should be removed?
	FULLY_INSIDE,  //whole cell is inside
	FULLY_OUTSIDE, //whole cell is outside
	PARTLY_INSIDE, //any point is inside
	PARTLY_OUTSIDE,//any point is outside
	CROSS,         //cell crosses/is tangent to the contour
	NO_CROSS,      //cell does not contact contour
};
struct TSubstractCells: public HMCallback::ExecutorBase{
	HMCB_SET_PROCNAME("Remove cells");
	HMCB_SET_DEFAULT_DURATION(100);

	GridData _run(const GridData& base, const Contour::Tree& cont, SubstractCellsAlgo algo);
};
extern HMCallback::FunctionWithCallback<TSubstractCells> SubstractCells;


struct OptInscribe{
	OptInscribe(double buffer_size,
		bool inside,
		int fillalgo=0,
		bool keep_cont=true,
		double angle0=0)
	:buffer_size(buffer_size), inside(inside),
			fillalgo(fillalgo), keep_cont(keep_cont), angle0(angle0){}

	double buffer_size; // shift
	bool inside;        // leave internal(=true) or external(=false) cont area
	int fillalgo;       // 0-triangle, 1-recombined, 99-no
	bool keep_cont;     // whether to use cont contour nodes(=true) or make repartition
	double angle0;      // significant angle for contour partition (only for keep_cont=false)
};

struct TInscribeGrid: public HMCallback::ExecutorBase{
	HMCB_SET_PROCNAME("Inscribe Grid");
	HMCB_SET_DEFAULT_DURATION(100);

	GridData _run(const GridData& base, const Contour::Tree& cont, OptInscribe opt);
};
extern HMCallback::FunctionWithCallback<TInscribeGrid> InscribeGrid;


struct OptInsertConstraints{
	OptInsertConstraints(double buffer_size,
		int fillalgo=0,
		bool keep_cont=true,
		double angle0=0)
	: buffer_size(buffer_size), fillalgo(fillalgo),
	  keep_cont(keep_cont), angle0(angle0){}

	double buffer_size;  // shift
	bool keep_cont;      //use line constraint segmentation(=true) or repart(=false);
	double angle0;      // significant angle for contour partition
	int fillalgo;       // 0-triangle, 1-recombined, 99-no
};
struct TInsertConstraints: public HMCallback::ExecutorBase{
	HMCB_SET_PROCNAME("Insert constraints");
	HMCB_SET_DEFAULT_DURATION(100);

	GridData _run(const GridData& base,
			const vector<EdgeData>& cont,                //Line constraint
			const vector<std::pair<Point, double>>& pnt, //Point + recommended step (-1 for auto step)
			OptInsertConstraints opt);
	GridData _run(const GridData& base,
			const vector<EdgeData>& cont,
			OptInsertConstraints opt);
	GridData _run(const GridData& base,
			const vector<std::pair<Point, double>>& pnt,
			OptInsertConstraints opt);
};
extern HMCallback::FunctionWithCallback<TInsertConstraints> InsertConstraints;


struct TSubstractArea: public HMCallback::ExecutorBase{
	HMCB_SET_PROCNAME("Domain exclusion");
	HMCB_SET_DEFAULT_DURATION(100);

	GridData _run(const GridData& g1, const Contour::Tree& area, bool is_inner);
};
extern HMCallback::FunctionWithCallback<TSubstractArea> SubstractArea;





}}}

#endif

