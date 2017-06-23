#include "inscribe_grid.hpp"
#include "infogrid.hpp"
#include "healgrid.hpp"

using namespace HM2D;
using namespace HM2D::Grid;

HMCallback::FunctionWithCallback<Algos::TSubstractCells> Algos::SubstractCells;
HMCallback::FunctionWithCallback<Algos::TInscribeGrid> Algos::InscribeGrid;

GridData Algos::TSubstractCells::_run(const GridData& base, const Contour::Tree& cont, SubstractCellsAlgo algo){
	auto invert_extraction = [&](const CellData& extr)->CellData{
		aa::constant_ids_pvec(base.vcells, 0);
		aa::constant_ids_pvec(extr, 1);
		CellData ret;
		std::copy_if(base.vcells.begin(), base.vcells.end(), std::back_inserter(ret),
				[](shared_ptr<Cell> c){ return c->id == 0; });
		return ret;
	};

	CellData goodcells;

	auto cb1 = callback->bottom_line_subrange(95);
	switch (algo){
	case SubstractCellsAlgo::FULLY_INSIDE:
	{
		CellData badcells = ExtractCells.WithCallback(cb1, base, cont, INSIDE);
		goodcells = invert_extraction(badcells);
		break;
	}
	case SubstractCellsAlgo::FULLY_OUTSIDE:
	{
		CellData badcells = ExtractCells.WithCallback(cb1, base, cont, OUTSIDE);
		goodcells = invert_extraction(badcells);
		break;
	}
	case SubstractCellsAlgo::PARTLY_INSIDE:
	{
		goodcells = ExtractCells.WithCallback(cb1, base, cont, OUTSIDE);
		break;
	}
	case SubstractCellsAlgo::PARTLY_OUTSIDE:
	{
		goodcells = ExtractCells.WithCallback(cb1, base, cont, INSIDE);
		break;
	}
	case SubstractCellsAlgo::CROSS:
	{
		CellData badcells = ExtractCells.WithCallback(cb1, base, cont, BOUND);
		goodcells = invert_extraction(badcells);
		break;
	}
	case SubstractCellsAlgo::NO_CROSS:
	{
		goodcells = ExtractCells.WithCallback(cb1, base, cont, BOUND);
		break;
	}
	}

	callback->step_after(5, "Assembling grid");
	GridData ret;
	DeepCopy(goodcells, ret.vcells, 2);
	RestoreFromCells(ret);

	return ret;
}

GridData Algos::TInscribeGrid::_run(const GridData& base, const Contour::Tree& cont,
		const OptInscribe& opt){

}
