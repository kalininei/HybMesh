#ifndef HYBMESH_DOMAPPING_HPP
#define HYBMESH_DOMAPPING_HPP

#include "gridmap.hpp"
#include "hmfem.hpp"
#include "mapped_contour.hpp"

namespace HMMap{namespace Impl{

struct DoMapping{
	DoMapping(const HMMap::Options& _opt, bool reversed):
		opt(_opt), mcol(reversed){}
	
	// ====== set input
	void set_grid(const HM2D::GridData& ig);
	void set_contour(const HM2D::EdgeData& ecol);
	void set_points(const vector<Point>& gridpnt, const vector<Point>& contpnt);

	// ====== main procedure
	//uses 100 units of cb
	HM2D::GridData run(HMCallback::Caller2& cb);
protected:
	// ====== main data
	HMMap::Options opt;
	HM2D::GridData inpgrid;
	HM2D::EdgeData contdata;
	vector<Point> gridpoints;
	vector<Point> contpoints;

	// ====== aux data
	//filled by prepare_mapped_contour
	HM2D::Contour::Tree mapped_outer;
	HM2D::Contour::Tree inpgrid_outer;
	//those are filled by prepare_grid and its subroutines
	shared_ptr<HM2D::GridData> g3;   //triangle grid
	HM2D::Contour::Tree g3outer;  //triangle grid borders
	shared_ptr<HMFem::LaplaceProblem> laplace; //laplace problem

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
	DirectMapping(Options opt, bool reversed): DoMapping(opt, reversed){}
private:
	void solve_uv_problems(vector<double>& u, vector<double>& v) override;
	//prepare_grid subroutines
	void build_grid3() override;
};

struct InverseMapping: public DoMapping{
	InverseMapping(Options opt, bool reversed): DoMapping(opt, reversed){}
private:
	void solve_uv_problems(vector<double>& u, vector<double>& v) override;
	//prepare_grid subroutines
	void build_grid3() override;
};

struct DoMapping2{
	DoMapping2(const HMMap::Options& _opt, bool reversed):
		opt(_opt), mcol(reversed){}
	
	// ====== set input
	void set_grid(const HM2D::GridData& ig);
	void set_contour(const HM2D::EdgeData& ecol);
	void set_points(const vector<Point>& gridpnt, const vector<Point>& contpnt);

	// ====== main procedure
	HM2D::GridData run();
private:
	// ====== main data
	HMMap::Options opt;
	HM2D::GridData inpgrid;
	HM2D::EdgeData contdata;
	vector<Point> gridpoints;
	vector<Point> contpoints;

	// ====== aux data
	//filled by prepare_mapped_contour
	HM2D::Contour::Tree mapped_outer;
	HM2D::Contour::Tree inpgrid_outer;
	//those are filled by prepare_grid and its subroutines
	shared_ptr<HM2D::GridData> g3;   //triangle grid within mapped_outer
	HM2D::Contour::Tree g3outer;  //triangle grid borders
	shared_ptr<HMFem::LaplaceProblem> laplace; //laplace problem

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
