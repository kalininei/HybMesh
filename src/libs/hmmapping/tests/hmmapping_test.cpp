#include <iostream>
#include "procgrid.h"
#include "hmmapping.hpp"
#include "hmtesting.hpp"
#include "circrect.hpp"
#include "vtk_export_grid2d.hpp"
using HMTesting::add_check;
using HMTesting::add_file_check;

void test01(){
	std::cout<<"01. Map grid from square"<<std::endl;
	GridGeom base = GGeom::Constructor::RectGrid01(10, 10);
	auto cont = HMCont2D::Constructor::ContourFromPoints({0,0, 1,0, 1.5, 1, 1, 1.3, 0.7, 1.0}, true);
	GridGeom mapped = HMGMap::MapGrid(base, cont,
		{Point(0, 0), Point(1, 0), Point(1,1), Point(0, 1)},
		{Point(0, 0), Point(1, 0), Point(1.5, 1), Point(0.7, 1.0)});
	add_check(base.n_cells() == mapped.n_cells(), "cells number");
	add_check(base.n_points() == mapped.n_points(), "nodes number");
	add_check([&](){
		double a1 = HMCont2D::Area(cont);
		double a2 = GGeom::Info::Area(mapped);
		if (a2>a1) return false;
		if (a1<0.9*a2) return false;
		return true;
	}(), "grids area");
}

void test02(){
	std::cout<<"02. Throws at invalid data"<<std::endl;
	GridGeom base = GGeom::Constructor::Ring(Point(1,1), 2, 1, 36, 6);

	auto c1 = HMCont2D::Constructor::Circle(4, 5, Point(0, 0));
	auto c2 = HMCont2D::Constructor::Circle(18, 2, Point(0, 0.1));
	HMCont2D::ECollection ecol;
	ecol.Unite(c1);
	ecol.Unite(c2);

	bool was_err1 = false;
	try{
		GridGeom ans1 = HMGMap::MapGrid(base, ecol,
			{Point(3, 1), Point(1, 3)},
			{Point(0, 5), Point(5, 0)}
		);
	} catch (HMGMap::MapException &e) {
		was_err1 = e.what() == (std::string)"Grid mapping exception: All grid boundaries should contain at least one base point";
	}
	add_check(was_err1, "grid boundary has no base points");

	bool was_err2 = false;
	try{
		GridGeom ans2 = HMGMap::MapGrid(base, ecol,
			{Point(3, 1), Point(1, 3), Point(2, 1)},
			{Point(0, 5), Point(5, 0), Point(0, -5)}
		);
	} catch (HMGMap::MapException &e) {
		was_err2 = e.what() == (std::string)"Grid mapping exception: Resulting grid is not valid";
	}
	add_check(was_err2, "points from same grid boundary are referenced to different mapped area");

	bool was_err3 = false;
	try{
		GridGeom ans3 = HMGMap::MapGrid(base, ecol,
			{Point(3, 1), Point(1, 3), Point(3, 3), Point(2, 1), Point(0, 1)},
			{Point(5, 0), Point(0, 5), Point(0, -5), Point(2, 0), Point(-2, 0)}
		);
	} catch (HMGMap::MapException &e) {
		was_err3 = e.what() == (std::string)"Grid mapping exception: Invalid order of points in mapped contour";
	}
	add_check(was_err3, "invalid points ordering");

	bool was_err4 = false;
	try{
		GridGeom ans4 = HMGMap::MapGrid(base, ecol,
			{Point(3, 1), Point(1, 3), Point(3, 3), Point(2, 1), Point(0, 1)},
			{Point(5, 0), Point(0, 5), Point(5, 5), Point(2, 0), Point(-2, 0)}
		);
	} catch (HMGMap::MapException &e) { was_err4 = true; }
	add_check(!was_err4, "valid data");
}

void test03(){
	std::cout<<"03. Circle vs rectangle"<<std::endl;
	GridGeom circgrid = GGeom::Constructor::Circle(Point(2, 3), 2, 16, 6, false);
	auto cont4 = HMCont2D::Constructor::ContourFromPoints({0,0, 3,0, 3,2, 0,2}, true);
	GridGeom rectgrid = GGeom::Constructor::RectGrid01(5, 7);
	auto contc = HMCont2D::Constructor::Circle(200, 16, Point(8, 9));

	GridGeom ans1 = HMGMap::MapGrid(circgrid, cont4,
			{Point(4, 3)},
			{Point(0, 0)}
	);
	add_check(ans1.n_points() == circgrid.n_points() &&
	          ans1.n_cells() == circgrid.n_cells(), "circle to rectangle, one base point");

	GridGeom ans2 = HMGMap::MapGrid(circgrid, cont4,
			{Point(4, 3), Point(2, 5)},
			{Point(0, 0), Point(3, 0)}
	);
	add_check(ans2.n_points() == circgrid.n_points() &&
	          ans2.n_cells() == circgrid.n_cells(), "circle to rectangle, two base points");

	GridGeom ans3 = HMGMap::MapGrid(circgrid, cont4,
			{Point(4, 3), Point(2, 5), Point(2, 1), Point(0, 3)},
			{Point(0, 0), Point(3, 0), Point(0, 2), Point(3,2)}
	);
	add_check(ans3.n_points() == circgrid.n_points() &&
	          ans3.n_cells() == circgrid.n_cells(), "circle to rectangle, four base points");

	GridGeom ans4 = HMGMap::MapGrid(rectgrid, contc,
	                {Point(0, 0)},
	                {Point(-8, 9)}
	);
	add_check(ans4.n_points() == rectgrid.n_points() &&
	          ans4.n_cells() == rectgrid.n_cells(), "rectangle to circle, one base point");
}

void test04(){
	std::cout<<"04. Snapping"<<std::endl;
	GridGeom rectgrid = GGeom::Constructor::RectGrid01(10, 10);
	auto mcont = HMCont2D::Constructor::ContourFromPoints(
		{0,0, 2,2, 1,3, -4,3}, true);
	HMGMap::Options opt;
	opt.fem_nrec = 5000;

	opt.snap = "ADD_VERTICES";
	GridGeom ans1 = HMGMap::MapGrid(rectgrid, mcont, {Point(0,0)}, {Point(0,0)}, opt);
	add_check(fabs(HMCont2D::Area(mcont) - GGeom::Info::Area(ans1))<1e-12 &&
	          rectgrid.n_points() + 3 == ans1.n_points(),
		"snapping by adding points");

	opt.snap = "SHIFT_VERTICES";
	GridGeom ans2 = HMGMap::MapGrid(rectgrid, mcont, {Point(0,0)}, {Point(0,0)}, opt);
	add_check(fabs(HMCont2D::Area(mcont) - GGeom::Info::Area(ans2))<1e-12 &&
	          rectgrid.n_points() == ans2.n_points(),
		"snapping by shifting points");
}

void test05(){
	std::cout<<"05. Circle with rectangle cells"<<std::endl;
	GridGeom g1 = HMGMap::Circ4Prototype(Point(0, 0), 1.0, 24, 1.0, 1.0);
	GGeom::Export::GridVTK(g1, "g1.vtk");
	auto sk = GGeom::Info::Skewness(g1);
	auto maxel = std::max_element(sk.begin(), sk.end());
	add_check(ISZERO(*maxel - 0.5), "skewness");
	add_file_check(8849316823715813599U, "g1.vtk", "no refinement, side = 1.0*rad");

	GridGeom g2 = HMGMap::Circ4Prototype(Point(3, 1), 2.0, 24, 1.0, 0.3);
	GGeom::Export::GridVTK(g2, "g1.vtk");
	add_file_check(7197682274436151092U, "g1.vtk", "with refinement, side = 1.0*rad");

	GridGeom g3 = HMGMap::Circ4Prototype(Point(-3, 1), 0.2, 80, 0.5, 1.0);
	GGeom::Export::GridVTK(g3, "g1.vtk");
	add_file_check(14945723505508438707U, "g1.vtk", "side = 0.5*rad");
}

int main(){
	test01();
	test02();
	test03();
	test04();
	test05();

	HMTesting::check_final_report();
	std::cout<<"DONE"<<std::endl;
	return 0;
}

