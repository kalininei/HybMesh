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
		std::vector<Point>{Point(0, 0), Point(1, 0), Point(1,1), Point(0, 1)},
		std::vector<Point>{Point(0, 0), Point(1, 0), Point(1.5, 1), Point(0.7, 1.0)},
		HMGMap::Options("inverse-laplace"));
	add_check(base.n_cells() == mapped.n_cells(), "cells number");
	add_check(base.n_points() == mapped.n_points(), "nodes number");
	add_check([&](){
		double a1 = HMCont2D::Area(cont);
		double a2 = GGeom::Info::Area(mapped);
		if (a2>a1) return false;
		if (a1<0.9*a2) return false;
		return true;
	}(), "grids area");
	GGeom::Export::GridVTK(mapped, "g1.vtk");
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
			std::vector<Point> {Point(3, 1), Point(1, 3)},
			std::vector<Point> {Point(0, 5), Point(5, 0)},
			HMGMap::Options("inverse-laplace")
		);
	} catch (HMGMap::MapException &e) {
		was_err1 = e.what() == (std::string)"Grid mapping exception: All grid boundaries should contain at least one base point";
	}
	add_check(was_err1, "grid boundary has no base points");

	bool was_err2 = false;
	try{
		GridGeom ans2 = HMGMap::MapGrid(base, ecol,
			std::vector<Point> {Point(3, 1), Point(1, 3), Point(2, 1)},
			std::vector<Point> {Point(0, 5), Point(5, 0), Point(0, -5)},
			HMGMap::Options("inverse-laplace")
		);
	} catch (HMGMap::MapException &e) {
		was_err2 = e.what() == (std::string)"Grid mapping exception: Contour-to-contour links are ambiguous";
	}
	add_check(was_err2, "points from same grid boundary are referenced to different mapped area");

	bool was_err3 = false;
	try{
		GridGeom ans3 = HMGMap::MapGrid(base, ecol,
			std::vector<Point> {Point(3, 1), Point(1, 3), Point(3, 3), Point(2, 1), Point(0, 1)},
			std::vector<Point> {Point(5, 0), Point(0, 5), Point(0, -5), Point(2, 0), Point(-2, 0)},
			HMGMap::Options("inverse-laplace")
		);
	} catch (HMGMap::MapException &e) {
		was_err3 = e.what() == (std::string)"Grid mapping exception: Invalid order of points in mapped contour";
	}
	add_check(was_err3, "invalid points ordering");

	try{
		HMGMap::Options opt("inverse-laplace");
		opt.fem_nrec = 1000;
		GridGeom ans4 = HMGMap::MapGrid(base, ecol,
			std::vector<Point> {Point(3, 1), Point(1, 3), Point(3, 3), Point(2, 1), Point(0, 1)},
			std::vector<Point> {Point(5, 0), Point(0, 5), Point(5, 5), Point(2, 0), Point(-2, 0)},
			opt
		);
		GGeom::Export::GridVTK(ans4, "g1.vtk");
		add_file_check(11959578928402032687U, "g1.vtk", "valid doubly connected data");
	} catch (HMGMap::MapException &e) {add_check(false, "valid doubly connected data");}
}

void test03(){
	std::cout<<"03. Circle vs rectangle"<<std::endl;
	GridGeom circgrid = GGeom::Constructor::Circle(Point(2, 3), 2, 16, 6, false);
	auto cont4 = HMCont2D::Constructor::ContourFromPoints({0,0, 3,0, 3,2, 0,2}, true);
	GridGeom rectgrid = GGeom::Constructor::RectGrid01(5, 7);
	auto contc = HMCont2D::Constructor::Circle(200, 16, Point(8, 9));

	GridGeom ans1 = HMGMap::MapGrid(circgrid, cont4,
			std::vector<Point> {Point(4, 3)},
			std::vector<Point> {Point(0, 0)},
			HMGMap::Options("inverse-laplace")
	);
	add_check(ans1.n_points() == circgrid.n_points() &&
	          ans1.n_cells() == circgrid.n_cells(), "circle to rectangle, one base point");

	GridGeom ans2 = HMGMap::MapGrid(circgrid, cont4,
			std::vector<Point> {Point(4, 3), Point(2, 5)},
			std::vector<Point> {Point(0, 0), Point(3, 0)},
			HMGMap::Options("inverse-laplace")
	);
	add_check(ans2.n_points() == circgrid.n_points() &&
	          ans2.n_cells() == circgrid.n_cells(), "circle to rectangle, two base points");

	GridGeom ans3 = HMGMap::MapGrid(circgrid, cont4,
			std::vector<Point> {Point(4, 3), Point(2, 5), Point(2, 1), Point(0, 3)},
			std::vector<Point> {Point(0, 0), Point(3, 0), Point(0, 2), Point(3,2)},
			HMGMap::Options("inverse-laplace")
	);
	add_check(ans3.n_points() == circgrid.n_points() &&
	          ans3.n_cells() == circgrid.n_cells(), "circle to rectangle, four base points");

	GridGeom ans4 = HMGMap::MapGrid(rectgrid, contc,
	                std::vector<Point> {Point(0, 0)},
	                std::vector<Point> {Point(-8, 9)},
			HMGMap::Options("inverse-laplace")
	);
	add_check(ans4.n_points() == rectgrid.n_points() &&
	          ans4.n_cells() == rectgrid.n_cells(), "rectangle to circle, one base point");
}

void test04(){
	std::cout<<"04. Snapping"<<std::endl;
	GridGeom rectgrid = GGeom::Constructor::RectGrid01(10, 10);
	auto mcont = HMCont2D::Constructor::ContourFromPoints(
		{0,0, 2,2, 1,3, -4,3}, true);
	HMGMap::Options opt("inverse-laplace");
	opt.fem_nrec = 5000;

	opt.snap = "ADD_VERTICES";
	GridGeom ans1 = HMGMap::MapGrid(rectgrid, mcont,
			std::vector<Point> {Point(0, 0)},
			std::vector<Point> {Point(0, 0)}, opt);
	add_check(fabs(HMCont2D::Area(mcont) - GGeom::Info::Area(ans1))<1e-12 &&
	          rectgrid.n_points() + 3 == ans1.n_points(),
		"snapping by adding points");
	GGeom::Export::GridVTK(ans1, "ans1.vtk");

	opt.snap = "SHIFT_VERTICES";
	GridGeom ans2 = HMGMap::MapGrid(rectgrid, mcont, vector<Point>{Point(0,0)}, vector<Point> {Point(0,0)}, opt);
	add_check(fabs(HMCont2D::Area(mcont) - GGeom::Info::Area(ans2))<1e-12 &&
	          rectgrid.n_points() == ans2.n_points(),
		"snapping by shifting points");
	GGeom::Export::GridVTK(ans2, "ans2.vtk");
}

void test05(){
	std::cout<<"05. Circle with rectangle cells"<<std::endl;
	GridGeom g1 = HMGMap::Circ4Prototype(Point(0, 0), 1.0, 24, 1.0, 1.0);
	GGeom::Export::GridVTK(g1, "g1.vtk");
	auto sk = GGeom::Info::Skewness(g1);
	auto maxel = std::max_element(sk.begin(), sk.end());
	add_check(ISZERO(*maxel - 0.5), "skewness");
	add_file_check(4015375048829588252U, "g1.vtk", "no refinement, side = 1.0*rad");

	GridGeom g2 = HMGMap::Circ4Prototype(Point(3, 1), 2.0, 24, 1.0, 0.3);
	GGeom::Export::GridVTK(g2, "g1.vtk");
	add_file_check(17078135237454258052U, "g1.vtk", "with refinement, side = 1.0*rad");

	GridGeom g3 = HMGMap::Circ4Prototype(Point(-3, 1), 0.2, 80, 0.5, 1.0);
	GGeom::Export::GridVTK(g3, "g1.vtk");
	add_file_check(14945723505508438707U, "g1.vtk", "side = 0.5*rad");
}

void test06(){
	std::cout<<"06. M-like target"<<std::endl;
	GridGeom g1 = GGeom::Constructor::RectGrid01(10, 10);
	auto cont = HMCont2D::Constructor::ContourFromPoints({0,0, 1,0, 1,1, 0.5, 0.2, 0,1}, true);

	auto ans = HMGMap::MapGrid.Silent(g1, cont,
			std::vector<Point> {Point(0,0), Point(1,0), Point(1,1), Point(0,1)},
			std::vector<Point> {Point(0,0), Point(1,0), Point(1,1), Point(0,1)},
			HMGMap::Options("inverse-laplace")
	);
	GGeom::Export::GridVTK(ans, "g1.vtk");
	add_file_check(5133065268170907878U, "g1.vtk", "m-like target, inverse");

	auto g1cont = GGeom::Info::Contour1(g1);
	auto ans2 = HMGMap::MapGrid.Silent(ans, g1cont,
			std::vector<Point> {Point(0,0), Point(1,0), Point(1,1), Point(0,1)},
			std::vector<Point> {Point(0,0), Point(1,0), Point(1,1), Point(0,1)},
			HMGMap::Options("direct-laplace")
	);
	GGeom::Export::GridVTK(ans2, "g1.vtk");
	add_file_check(4588482714625395358U, "g1.vtk", "m-like base, direct");
}

void test07(){
	std::cout<<"07. Orthonal grid in curvilinear quadrangle"<<std::endl;
	{
		//initial contours
		HMCont2D::PCollection pcol;
		auto left1 = HMCont2D::Constructor::PerturbedContour(Point(0, 0), Point(-0.11, 1), 1000,
				[](double x){ return 0.2*sin(2*M_PI*x); });
		auto bot1 = HMCont2D::Constructor::PerturbedContour(Point(0, 0), Point(3, 0), 1000,
				[](double x){ return 0.13*sin(8*M_PI*x); });
		auto right1 = HMCont2D::Constructor::PerturbedContour(Point(3, 0), Point(3.06, 1), 1000,
				[](double x){ return 0.07*sin(4*M_PI*x); });
		auto top1 = HMCont2D::Constructor::PerturbedContour(Point(-0.11, 1), Point(3.06, 1), 1000,
				[](double x){ return 0.1*sin(6*M_PI*x); });
		//partition
		auto left = HMCont2D::Algos::Partition(0.1, left1, pcol);
		auto bot = HMCont2D::Algos::Partition(0.1, bot1, pcol);
		auto right = HMCont2D::Algos::Partition(0.1, right1, pcol);
		auto top = HMCont2D::Algos::Partition(0.1, top1, pcol);
		//build grid
		HMCont2D::SaveVtk(left, "left.vtk");
		HMCont2D::SaveVtk(bot, "bot.vtk");
		HMCont2D::SaveVtk(right, "right.vtk");
		HMCont2D::SaveVtk(top1, "top.vtk");
		GridGeom ans = HMGMap::OrthogonalRectGrid(left, bot, right1, top1);
		GGeom::Export::GridVTK(ans, "g1.vtk");
	}
}


int main(){
	test01();
	test02();
	test03();
	test04();
	test05();
	test06();
	test07();

	HMTesting::check_final_report();
	std::cout<<"DONE"<<std::endl;
	return 0;
}

