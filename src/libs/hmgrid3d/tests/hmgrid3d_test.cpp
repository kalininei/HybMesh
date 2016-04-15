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

	add_check(g1.n_vertices() == 64 && g1.n_cells() == 27 &&
			g1.n_edges() == 144 && g1.n_faces() == 108,
			"cuboid primitives number"); 
	add_file_check(8514389607504677061U, "g1.vtk", "grid");
	add_file_check(12288405335072816377U, "c1.vtk", "boundary");
}

void test02(){
	std::cout<<"2. parallel sweep"<<std::endl;
	{
		auto g2d = GGeom::Constructor::RectGrid01(11, 7);
		auto g3d = HMGrid3D::Constructor::SweepGrid2D(g2d, {0.3, 0.4, 0.8});
		HMGrid3D::Export::BoundaryVTK(g3d, "c1.vtk");
		HMGrid3D::Export::GridVTK(g3d, "g1.vtk");
		add_check(g3d.n_vertices() == 288 && g3d.n_edges() == 708 &&
				g3d.n_faces() == 575 && g3d.n_cells() == 154,
				"rectangular grid sweep");
	}
	{
		auto g2d = GGeom::Constructor::Ring(Point(0, 0), 4, 2, 12, 4);
		auto g3d = HMGrid3D::Constructor::SweepGrid2D(g2d, {0, 0.1, 0.4, 0.7});
		HMGrid3D::Export::BoundaryVTK(g3d, "c1.vtk");
		HMGrid3D::Export::GridVTK(g3d, "g1.vtk");
		add_check(g3d.n_vertices() == 240 && g3d.n_cells() == 144,
				"ring grid sweep");
	}
	{
		auto g2d = GGeom::Constructor::Circle(Point(1, 0), 4, 24, 10, true);
		auto g3d = HMGrid3D::Constructor::SweepGrid2D(g2d, {0, 1, 2, 3});
		HMGrid3D::Export::BoundaryVTK(g3d, "c1.vtk");
		HMGrid3D::Export::GridVTK(g3d, "g1.vtk");
		add_check(g3d.n_vertices() == 964 && g3d.n_cells() == 720,
				"circle grid sweep");
	}
}

void test03(){
	std::cout<<"3. Fluent export"<<std::endl;
	{
		auto g1 = HMGrid3D::Constructor::Cuboid(HMGrid3D::Vertex(0, 0, 0), 1, 1, 1, 2, 2, 1);
		HMGrid3D::Export::GridMSH.Silent(g1, "g1.msh");
		add_file_check(11216638150753088583U, "g1.msh", "simple cuboid");
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
		add_file_check(1976956498869484852U, "g1.msh", "cuboid from sweep with custom boundaries");
	}
	{
		auto g1 = GGeom::Constructor::Circle(Point(0, 0), 1, 4, 2, true);
		auto g2 = GGeom::Constructor::ExtractCells(g1, {0, 4});
		auto g3d = HMGrid3D::Constructor::SweepGrid2D(g2, {0, 0.1});
		HMGrid3D::Export::GridMSH(g3d, "g1.msh");
		add_file_check(17687514312539715596U, "g1.msh", "mixed hex/wedge cells");
	}
	{
		auto g1 = GGeom::Constructor::Circle(Point(0, 0), 1, 5, 2, false);
		auto g2 = GGeom::Constructor::ExtractCells(g1, {5});
		auto g3d = HMGrid3D::Constructor::SweepGrid2D(g2, {0, 0.1});
		HMGrid3D::Export::GridMSH(g3d, "g1.msh");
		add_file_check(3921094917079304U, "g1.msh", "single pentagon prism cell");
	}
	{
		auto g1 = GGeom::Constructor::RectGrid01(20, 30);
		auto g2 = GGeom::Constructor::Circle(Point(0.721, 0.682), 0.465, 24, 10, false);
		auto g3 = GridGeom::cross_grids(&g1, &g2, 0.0, 0, false, false, 0);
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
		add_file_check(1949417467069346897U, "g1.msh", "mesh with polyhedra cells");
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
		add_file_check(17371627965264670427U, "_o1.msh", "simple 2x2x2");

		pd.data[0].reversed = false;
		HMGrid3D::Export::GridMSH(g3d, "_o2.msh", pd); 
		add_file_check(9374713153386878447U, "_o2.msh", "2x2x2 without reverse");

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
		auto g2d = GridGeom::cross_grids(&g2d1, &g2d2, 0.1, 0, false, true, 0);
		vector<double> zvec;
		for (int i=0; i<100; i+=10)  zvec.push_back(3 + (double)i/99);
		auto g3d = HMGrid3D::Constructor::SweepGrid2D(*g2d, zvec);
		HMGrid3D::Face::SetBoundaryTypes(g3d.allfaces(), [](HMGrid3D::Vertex v, int bt){
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
		add_file_check(6712434711251990504U, "g2.msh", "multiple periodic");

		delete g2d;
	}
};

void test05(){
	std::cout<<"5. Tecplot export"<<std::endl;
	{
		auto g2d = GGeom::Constructor::RectGrid01(1, 1);
		auto g3d = HMGrid3D::Constructor::SweepGrid2D(g2d, {0.0, 0.5});
		HMGrid3D::Export::GridTecplot.Silent(g3d, "g1.dat");
		add_file_check(1996510299573747148U, "g1.dat", "single cell grid");
	}
	{
		auto g2d = GGeom::Constructor::Circle(Point(0, 0), 10, 30, 10, false);
		auto g3d = HMGrid3D::Constructor::SweepGrid2D(g2d, {1.0, 1.2, 1.4, 1.6, 1.7, 1.8, 1.9, 2.0});
		HMGrid3D::Export::GridTecplot.Silent(g3d, "g1.dat");
		add_file_check(15174883252499084243U, "g1.dat", "polyhedral grid");
		HMGrid3D::Export::BoundaryTecplot.Silent(g3d, "g1.dat");
		add_file_check(4351852141727252628U, "g1.dat", "polyhedral boundary");
	}
}
int main(){
	test01();
	test02();
	test03();
	test04();
	test05();
	
	check_final_report();
	std::cout<<"DONE"<<std::endl;
}
