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
	OptInscribe(double buffer_size, bool inside,
		int fillalgo=0,
		bool keep_cont=true,
		double angle0=0)
	:buffer_size(buffer_size), inside(inside),
		fillalgo(fillalgo), keep_cont(keep_cont), angle0(angle0){};

	double buffer_size; // shift
	bool inside;        // leave internal(=true) or external(=false) cont area
	int fillalgo;       // 0-triangle, 1-recombined, 99-no
	bool keep_cont;     // whether to use contour nodes(=true) or make contour partition
	double angle0;      // significant angle for contour partition (only for keep_cont=true)
};

struct TInscribeGrid: public HMCallback::ExecutorBase{
	HMCB_SET_PROCNAME("Inscribe Grid");
	HMCB_SET_DEFAULT_DURATION(100);

	GridData _run(const GridData& base, const Contour::Tree& cont, OptInscribe opt);
};
extern HMCallback::FunctionWithCallback<TInscribeGrid> InscribeGrid;


}}}

#endif

