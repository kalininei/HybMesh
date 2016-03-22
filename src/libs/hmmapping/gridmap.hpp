#ifndef HYBMESH_GRIDMAP_HPP
#define HYBMESH_GRIDMAP_HPP

#include "procgrid.h"
#include "hmmapping.hpp"
#include "femgrid43.hpp"
#include "hmfem.hpp"
#include "mapped_contour.hpp"

namespace HMGMap{
namespace Impl{


struct DoMapping{
	DoMapping(const Options& _opt): opt(_opt), inpgrid(GGeom::Constructor::EmptyGrid()){}
	
	// ====== set input
	void set_grid(const GridGeom& ig);
	void set_contour(const HMCont2D::ECollection& ecol);
	void set_points(const vector<Point>& gridpnt, const vector<Point>& contpnt);

	// ====== main procedure
	GridGeom run();
private:

	// ====== main data
	Options opt;
	GridGeom inpgrid;
	HMCont2D::ECollection contdata;
	vector<Point> gridpoints;
	vector<Point> contpoints;

	// ====== aux data
	//filled by prepare_mapped_contour
	HMCont2D::ContourTree mapped_outer;
	//those are filled by prepare_grid and its subroutines
	shared_ptr<HMFem::Grid43> g3;   //triangle grid
	HMCont2D::ContourTree g3outer;  //triangle grid borders
	shared_ptr<HMFem::Mat> laplas_mat;  //preassembled laplas operator
	MappedContourCollection mcol;   //contour mapper: from g3 contour to contdata contour

	//aux procedures
	void prepare_mapped_contour();
	void prepare_grid();
	vector<double> solve_u_problem();
	vector<double> solve_v_problem();
	//prepare_grid subroutines
	void build_grid3();
	void build_mcc();

};



}}


#endif
