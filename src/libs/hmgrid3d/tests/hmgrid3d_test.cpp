#include "hmgrid3d.hpp"
#include "procgrid.h"
#include <fstream>
#include "debug_grid3d.hpp"
#include "debug_grid2d.h"
#include "hmtesting.hpp"
#include "hmtimer.hpp"
using namespace HMTesting;

void test01(){
	std::cout<<"1. export cuboid to vtk"<<std::endl;
	auto g1 = HMGrid3D::Constructor::Cuboid({0, 0, 0}, 1, 2, 5, 3, 3, 3);
	HMGrid3D::Export::AllVTK.Silent(g1, "g1.vtk", "c1.vtk");
	add_check(g1.n_vert == 64 && g1.n_cells == 27 &&
			g1.n_edges == 144 && g1.n_faces == 108,
			"cuboid primitives number"); 
	add_file_check(15732503264642486832U, "g1.vtk", "grid");
	add_file_check(12574868808136614456U, "c1.vtk", "boundary");
}

void test02(){
	std::cout<<"2. parallel sweep"<<std::endl;
	{
		auto g2d = GGeom::Constructor::RectGrid01(11, 7);
		auto g3d = HMGrid3D::Constructor::SweepGrid2D(g2d, {0.3, 0.4, 0.8});
		HMGrid3D::Export::BoundaryVTK(g3d, "c1.vtk");
		HMGrid3D::Export::GridVTK(g3d, "g1.vtk");
		add_check(g3d.n_vert == 288 && g3d.n_edges == 708 &&
				g3d.n_faces == 575 && g3d.n_cells == 154,
				"rectangular grid sweep");
	}
	{
		auto g2d = GGeom::Constructor::Ring(Point(0, 0), 4, 2, 12, 4);
		auto g3d = HMGrid3D::Constructor::SweepGrid2D(g2d, {0, 0.1, 0.4, 0.7});
		HMGrid3D::Export::BoundaryVTK(g3d, "c1.vtk");
		HMGrid3D::Export::GridVTK(g3d, "g1.vtk");
		add_check(g3d.n_vert == 240 && g3d.n_cells == 144,
				"ring grid sweep");
	}
	{
		auto g2d = GGeom::Constructor::Circle(Point(1, 0), 4, 24, 10, true);
		auto g3d = HMGrid3D::Constructor::SweepGrid2D(g2d, {0, 1, 2, 3});
		HMGrid3D::Export::BoundaryVTK(g3d, "c1.vtk");
		HMGrid3D::Export::GridVTK(g3d, "g1.vtk");
		add_check(g3d.n_vert == 964 && g3d.n_cells == 720,
				"circle grid sweep");
	}
}

void test03(){
	std::cout<<"3. Fluent export"<<std::endl;
	{
		auto g1 = HMGrid3D::Constructor::Cuboid(HMGrid3D::Vertex(0, 0, 0), 1, 1, 1, 2, 2, 1);
		HMGrid3D::Export::GridMSH.Silent(g1, "g1.msh");
		add_file_check(4712463904373753464U, "g1.msh", "simple cuboid");
	}
	{
		auto g2d = GGeom::Constructor::RectGrid01(6, 3);
		auto g3d = HMGrid3D::Constructor::SweepGrid2D(g2d, {0, 0.1, 0.2, 0.5},
				[](int i){ return 1; },
				[](int i){ return 2; },
				[](int i){ return i+3; });
		HMGrid3D::Export::GridMSH(g3d, "g1.msh",
				[](int i)->std::string{
					switch (i){
						case 1: return std::string("bottom");
						case 2: return std::string("top");
						default: return std::string("side") + std::to_string(i);
					};
				});
		add_file_check(15070544857795696967U, "g1.msh", "cuboid from sweep with custom boundaries");
	}
	{
		auto g1 = GGeom::Constructor::Circle(Point(0, 0), 1, 4, 2, true);
		auto g2 = GGeom::Constructor::ExtractCells(g1, {0, 4});
		auto g3d = HMGrid3D::Constructor::SweepGrid2D(g2, {0, 0.1});
		HMGrid3D::Export::GridMSH(g3d, "g1.msh");
		add_file_check(14999299632248408018U, "g1.msh", "mixed hex/wedge cells");
	}
	{
		auto g1 = GGeom::Constructor::Circle(Point(0, 0), 1, 5, 2, false);
		auto g2 = GGeom::Constructor::ExtractCells(g1, {5});
		auto g3d = HMGrid3D::Constructor::SweepGrid2D(g2, {0, 0.1});
		HMGrid3D::Export::GridMSH(g3d, "g1.msh");
		add_file_check(17235633717941475062U, "g1.msh", "single pentagon prism cell");
	}
	{
		auto g1 = GGeom::Constructor::RectGrid01(20, 30);
		auto g2 = GGeom::Constructor::Circle(Point(0.721, 0.682), 0.465, 24, 10, false);
		auto g3 = GridGeom::cross_grids(&g1, &g2, 0.0, 0, false, false, 0, 0);
		auto _g3ed = g3->get_edges();
		vector<::Edge> g3ed(_g3ed.begin(), _g3ed.end());
		auto g3d = HMGrid3D::Constructor::SweepGrid2D(*g3, {0, 0.1, 0.2, 0.3, 0.5},
				[](int){return 1;}, [](int){return 2;},
				[&g3, &g3ed](int e)->int{
					Point pc = (*g3->get_point(g3ed[e].p1) + *g3->get_point(g3ed[e].p2))/2.0;
					return (pc.x<=1+1e-12 && pc.y<=1+1e-12) ? 3 : 4;
				});
		HMGrid3D::Export::GridMSH(g3d, "g1.msh",
				[](int i)->std::string{
					switch (i){
						case 1: return "bottom";
						case 2: return "top";
						case 3: return "square";
						case 4: return "circle";
						default: return "unknown";
					}
				});
		add_file_check(9690145692750233246U, "g1.msh", "mesh with polyhedra cells");
		delete g3;
	}
}

void test04(){
	std::cout<<"4. Fluent export with periodic surfaces"<<std::endl;
	{
		auto g2d = GGeom::Constructor::RectGrid({0.0, 0.1, 1.0}, {0.0, 0.3, 1.0});
		auto g3d = HMGrid3D::Constructor::SweepGrid2D(g2d, {0.0, 0.2, 1.0});
		HMGrid3D::Export::PeriodicData pd;
		pd.add_condition(1, 2, HMGrid3D::Vertex(0, 0, 0), HMGrid3D::Vertex(0, 0, 1), true);
		HMGrid3D::Export::GridMSH(g3d, "_o1.msh", pd); 
		add_file_check(9831086485199871919U, "_o1.msh", "simple 2x2x2");

		pd.data[0].reversed = false;
		HMGrid3D::Export::GridMSH(g3d, "_o2.msh", pd); 
		add_file_check(7770892506734403513U, "_o2.msh", "2x2x2 without reverse");

		pd.data[0].reversed = true;
		pd.data[0].v = HMGrid3D::Vertex(0.1, 0, 0);
		bool res3 = false;
		try{
			HMGrid3D::Export::GridMSH(g3d, "g3.msh", pd); 
		} catch(std::runtime_error& e){
			res3 = true;
		}
		add_check(res3, "fail at invalid point match");
	}
	{
		auto g2d1 = GGeom::Constructor::RectGrid(Point(0, 0), Point(10, 1), 100, 10);
		auto g2d2 = GGeom::Constructor::Ring(Point(3, 0.5), 0.3, 0.1, 20, 4);
		auto g2d = GridGeom::cross_grids(&g2d1, &g2d2, 0.1, 0, false, true, 0, 0);
		vector<double> zvec;
		for (int i=0; i<100; i+=10)  zvec.push_back(3 + (double)i/99);
		auto g3d = HMGrid3D::Constructor::SweepGrid2D(*g2d, zvec);
		g3d.set_btype([](HMGrid3D::Vertex v, int bt){
					if (bt == 3){
						if (v.x<=geps) return 3;
						if (v.x>=10-geps) return 4;
						return 5;
					}
					return bt;
				});

		HMGrid3D::Export::PeriodicData pd;

		pd.add_condition(1, 2, HMGrid3D::Vertex(0, 0, 3), HMGrid3D::Vertex(0, 0, 4), true);
		pd.add_condition(3, 4, HMGrid3D::Vertex(0, 0, 3), HMGrid3D::Vertex(10, 0, 3), true);
		HMGrid3D::Export::GridMSH.Silent(g3d, "g2.msh", pd);
		add_file_check(16533597872948878492U, "g2.msh", "multiple periodic");
		//#####
		throw;

		delete g2d;
	}
};

void test05(){
	std::cout<<"5. Tecplot export"<<std::endl;
	{
		auto g2d = GGeom::Constructor::RectGrid01(1, 1);
		auto g3d = HMGrid3D::Constructor::SweepGrid2D(g2d, {0.0, 0.5});
		HMGrid3D::Export::GridTecplot.Silent(g3d, "g1.dat");
		add_file_check(1831833575709478659U, "g1.dat", "single cell grid");
	}
	{
		auto g2d = GGeom::Constructor::Circle(Point(0, 0), 10, 30, 10, false);
		auto g3d = HMGrid3D::Constructor::SweepGrid2D(g2d, {1.0, 1.2, 1.4, 1.6, 1.7, 1.8, 1.9, 2.0});
		HMGrid3D::Export::GridTecplot.Silent(g3d, "g1.dat");
		add_file_check(17626851046985520587U, "g1.dat", "polyhedral grid");
		HMGrid3D::Export::BoundaryTecplot.Silent(g3d, "g1.dat");
		add_file_check(8291026423155100327U, "g1.dat", "polyhedral boundary");
	}
}

void test06(){
	using HMGrid3D::Constructor::RevolveGrid2D;
	std::cout<<"6. Solid of revolution"<<std::endl;
	auto g2d = GGeom::Constructor::RectGrid(Point(1,0), Point(2,1), 1, 1);
	auto bc0 = [](int){return 0;};
	{
		auto g3d = RevolveGrid2D(g2d, {0, 90}, Point(0, 0), Point(0, 1), true, bc0, bc0, bc0);
		HMGrid3D::Export::GridTecplot.Silent(g3d, "g1.dat");
		add_file_check(16088294825526263046U, "g1.dat", "single cell, distant, incomplete");
	}
	{
		auto g3d = RevolveGrid2D(g2d, {0, 90, 180, 270, 360}, Point(0, 0), Point(0, 1), true,
				bc0, bc0, bc0);
		HMGrid3D::Export::GridTecplot.Silent(g3d, "g1.dat");
		add_file_check(8732440237994901672U, "g1.dat", "single cell, distant, complete");
	}
	{
		auto g3d = RevolveGrid2D(g2d, {0, 90, 100}, Point(0, 0), Point(0, 1), true,
				[](int i){ return i; }, [](int){return 10;}, [](int){return 20;});
		HMGrid3D::Export::GridTecplot.Silent(g3d, "g1.dat");
		add_file_check(3859847262675033285U, "g1.dat", "single cell, distant, incomplete, with bc");
	}
	{
		auto h2d = GGeom::Constructor::RectGrid(Point(0, 0), Point(2, 1), 2, 1);
		auto g3d = RevolveGrid2D(h2d, {0, 90}, Point(0, 0), Point(0, 1), true);
		HMGrid3D::Export::GridTecplot.Silent(g3d, "g1.dat");
		add_file_check(8233442907809870919U, "g1.dat", "with contact, incomplete");
	}
	{
		auto h2d = GGeom::Constructor::RectGrid(Point(0, 0), Point(2, 1), 4, 3);
		auto g3d = RevolveGrid2D(h2d, {0, 90, 110, 180, 250, 330, 360}, Point(0, 0), Point(0, 1), true);
		HMGrid3D::Export::GridTecplot.Silent(g3d, "g1.dat");
		add_file_check(5490115627065179709U, "g1.dat", "with contact, complete");
	}
	{
		auto g1 = GGeom::Constructor::RectGrid(Point(0, 0), Point(10, 10), 10, 10);
		auto g2 = GGeom::Constructor::RectGrid(Point(0, 5), Point(10, 6), 5, 1);
		auto g3 = GridGeom::cross_grids(&g1, &g2, 0.0, 0, 0, 0, 0, 0);
		auto g3d = RevolveGrid2D(*g3, {0, 10, 20, 30}, Point(0, 0), Point(0, 1), true);
		HMGrid3D::Export::GridTecplot.Silent(g3d, "g1.dat");
		add_file_check(12980710001405184230U, "g1.dat", "hanging nodes near axis to tecplot");
		HMGrid3D::Export::GridMSH.Silent(g3d, "g1.msh");
		add_file_check(7731290871669899787U, "g1.msh", "hanging nodes near axis to fluent");
		delete g3;
	}
}

void test07(){
	using HMGrid3D::Constructor::RevolveGrid2D;
	std::cout<<"7. Solid of revolution, merging centeral cells"<<std::endl;
	{
		auto g2d = GGeom::Constructor::RectGrid(Point(1,0), Point(2,1), 1, 1);
		auto g3d = RevolveGrid2D(g2d, {0, 45, 90}, Point(1, 0), Point(1, 1), false);
		HMGrid3D::Export::GridTecplot.Silent(g3d, "g1.dat");
		add_file_check(13398422286724743124U, "g1.dat", "single cell, without center trian, incomplete");
	}
	{
		auto g2d = GGeom::Constructor::RectGrid(Point(1,0), Point(2,1), 1, 1);
		auto g3d = RevolveGrid2D(g2d, {20, 45, 90, 160, 270, 300, 380}, Point(1, 0), Point(1, 1), false);
		HMGrid3D::Export::GridTecplot.Silent(g3d, "g1.dat");
		add_file_check(6994418583934313116U, "g1.dat", "single cell, without center trian, complete");
	}
	{
		auto h2d = GGeom::Constructor::RectGrid(Point(0, 0), Point(2, 1), 2, 1);
		auto g3d = RevolveGrid2D(h2d, {0, 10, 20, 30, 40, 50}, Point(0, 0), Point(0, 1), false);
		HMGrid3D::Export::GridTecplot.Silent(g3d, "g1.dat");
		add_file_check(11881236001573517783U, "g1.dat", "multiple cells, with trian, complete");
	}
	{
		auto g1 = GGeom::Constructor::EmptyGrid();
		GGeom::Modify::AddCell(g1, {Point(0,0), Point(1,0), Point(0, 1)});
		auto g2 = RevolveGrid2D(g1, {0, 45, 90}, Point(0,0), Point(0,1), false);
		HMGrid3D::Export::GridTecplot.Silent(g2, "g1.dat");
		add_file_check(10167032458429566145U, "g1.dat", "no tri with single axis triangle");
	}
	{
		auto g1 = GGeom::Constructor::EmptyGrid();
		GGeom::Modify::AddCell(g1, {Point(1,0), Point(1,1), Point(0, 1)});
		auto g2 = RevolveGrid2D(g1, {0, 45, 90}, Point(0,0), Point(0,1), false);
		HMGrid3D::Export::GridTecplot.Silent(g2, "g1.dat");
		add_file_check(11550191908304285294U, "g1.dat", "no tri, single off axis triangle");
	}
	{
		auto g1 = GGeom::Constructor::EmptyGrid();
		GGeom::Modify::AddCell(g1, {Point(0,0), Point(1,0), Point(0, 1)});
		GGeom::Modify::AddCell(g1, {Point(1,0), Point(1,1), Point(0, 1)});
		GGeom::Modify::AddCell(g1, {Point(0,0), Point(0,-2), Point(1, -2), Point(1, 0)});
		GGeom::Repair::Heal(g1);
		auto g2 = RevolveGrid2D(g1, {0, 45, 90}, Point(0,0), Point(0,1), false);
		HMGrid3D::Export::GridTecplot.Silent(g2, "g1.dat");
		add_file_check(12664340621499564857U, "g1.dat", "no tri, complex connections, incomplete");
		auto g3 = RevolveGrid2D(g1, {0, 90, 180, 270, 360}, Point(0,0), Point(0,1), false);
		HMGrid3D::Export::GridTecplot.Silent(g3, "g1.dat");
		add_file_check(2848618625331037303U, "g1.dat", "no tri, complex connections, complete");
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
	
	check_final_report();
	std::cout<<"DONE"<<std::endl;
}
