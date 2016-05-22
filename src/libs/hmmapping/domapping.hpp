#ifndef HYBMESH_DOMAPPING_HPP
#define HYBMESH_DOMAPPING_HPP

#include "hybmesh_contours2d.hpp"
#include "gridmap.hpp"
#include "hmfem.hpp"
#include "mapped_contour.hpp"

namespace HMMap{namespace Impl{

struct DoMapping{
	DoMapping(const HMMap::Options& _opt): opt(_opt), inpgrid(GGeom::Constructor::EmptyGrid()){}
	
	// ====== set input
	void set_grid(const GridGeom& ig);
	void set_contour(const HMCont2D::ECollection& ecol);
	void set_points(const vector<Point>& gridpnt, const vector<Point>& contpnt);

	// ====== main procedure
	//uses 100 units of cb
	GridGeom run(HMCallback::Caller2& cb);
protected:
	// ====== main data
	HMMap::Options opt;
	GridGeom inpgrid;
	HMCont2D::ECollection contdata;
	vector<Point> gridpoints;
	vector<Point> contpoints;

	// ====== aux data
	//filled by prepare_mapped_contour
	HMCont2D::ContourTree mapped_outer;
	HMCont2D::ContourTree inpgrid_outer;
	//those are filled by prepare_grid and its subroutines
	shared_ptr<HMFem::Grid43> g3;   //triangle grid
	HMCont2D::ContourTree g3outer;  //triangle grid borders
	shared_ptr<HMFem::LaplasProblem> laplace; //laplace problem

	MappedContourCollection mcol;   //contour mapper: from g3 contour to contdata contour

	//aux procedures
	void prepare_mapped_contour();
	void prepare_grid();
	virtual void solve_uv_problems(vector<double>& u, vector<double>& v) = 0;
	//prepare_grid subroutines
	virtual void build_grid3() = 0;
	void build_mcc();
};

struct DirectMapping: public DoMapping{
	DirectMapping(Options opt): DoMapping(opt){}
private:
	void solve_uv_problems(vector<double>& u, vector<double>& v) override;
	//prepare_grid subroutines
	void build_grid3() override;
};

struct InverseMapping: public DoMapping{
	InverseMapping(Options opt): DoMapping(opt){}
private:
	void solve_uv_problems(vector<double>& u, vector<double>& v) override;
	//prepare_grid subroutines
	void build_grid3() override;
};

struct DoMapping2{
	DoMapping2(const HMMap::Options& _opt): opt(_opt), inpgrid(GGeom::Constructor::EmptyGrid()){}
	
	// ====== set input
	void set_grid(const GridGeom& ig);
	void set_contour(const HMCont2D::ECollection& ecol);
	void set_points(const vector<Point>& gridpnt, const vector<Point>& contpnt);

	// ====== main procedure
	GridGeom run();
private:
	// ====== main data
	HMMap::Options opt;
	GridGeom inpgrid;
	HMCont2D::ECollection contdata;
	vector<Point> gridpoints;
	vector<Point> contpoints;

	// ====== aux data
	//filled by prepare_mapped_contour
	HMCont2D::ContourTree mapped_outer;
	HMCont2D::ContourTree inpgrid_outer;
	//those are filled by prepare_grid and its subroutines
	shared_ptr<HMFem::Grid43> g3;   //triangle grid within mapped_outer
	HMCont2D::ContourTree g3outer;  //triangle grid borders
	shared_ptr<HMFem::LaplasProblem> laplace; //laplace problem

	MappedContourCollection mcol;   //contour mapper: from g3 contour to inpgrid contour

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
