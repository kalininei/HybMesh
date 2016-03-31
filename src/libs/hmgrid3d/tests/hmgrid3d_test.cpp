#include "hmgrid3d.hpp"
#include "procgrid.h"
#include "fileproc.h"
#include <fstream>
#include "debug_grid3d.hpp"

int FAILED_CHECKS = 0;

void add_check(bool ex, std::string info){
	if (info.size()==0){
		std::cout<<"\tunknown check: ";
	} else{
		std::cout<<"\t"<<info;
	}
	if (ex){
		std::cout<<": True"<<std::endl;
	} else {
		++FAILED_CHECKS;
		std::cout<<": False <<<<<<<<<<<<<<<<<<<"<<std::endl;
	}
};

void add_file_check(long hash, std::string fn, std::string info){
	std::ifstream t(fn);
	std::stringstream buffer;
	buffer << t.rdbuf();
	std::hash<std::string> str_hash;
	long h1 = str_hash(buffer.str());
	int fc = FAILED_CHECKS;
	add_check(h1 == hash, info);
	if (fc != FAILED_CHECKS){
		std::cout<<"\t\tcurrent file hash = "<<h1<<std::endl;
		std::cout<<"\t\tinput file hash   = "<<hash<<std::endl;
	}
}

void test01(){
	std::cout<<"1. export import cuboid"<<std::endl;
	auto g1 = HMGrid3D::Constructor::Cuboid({0, 0, 0}, 1, 2, 5, 3, 3, 3);
	HMGrid3D::Export::BoundaryVTK(g1, "c1.vtk");
	HMGrid3D::Export::GridVTK(g1, "g1.vtk");
	add_check(g1.n_vertices() == 64 && g1.n_cells() == 27 &&
			g1.n_edges() == 144 && g1.n_faces() == 108,
			"cuboid primitives number"); 
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
		HMGrid3D::Export::GridMSH(g1, "g1.msh");
		add_file_check(-8255129009601487661, "g1.msh", "simple cuboid");
		 
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
		add_file_check(-7221922955056885602, "g1.msh", "cuboid from sweep with custom boundaries");
	}
	{
		auto g1 = GGeom::Constructor::Circle(Point(0, 0), 1, 4, 2, true);
		auto g2 = GGeom::Constructor::ExtractCells(g1, {0, 4});
		auto g3d = HMGrid3D::Constructor::SweepGrid2D(g2, {0, 0.1});
		HMGrid3D::Export::GridMSH(g3d, "g1.msh");
		add_file_check(-1907392421441205334, "g1.msh", "mixed hex/wedge cells");
	}
	{
		auto g1 = GGeom::Constructor::Circle(Point(0, 0), 1, 5, 2, false);
		auto g2 = GGeom::Constructor::ExtractCells(g1, {5});
		auto g3d = HMGrid3D::Constructor::SweepGrid2D(g2, {0, 0.1});
		HMGrid3D::Export::GridMSH(g3d, "g1.msh");
		add_file_check(-3354809844770748450, "g1.msh", "single pentagon prism cell");
	}
	{
		auto g1 = GGeom::Constructor::RectGrid01(20, 30);
		auto g2 = GGeom::Constructor::Circle(Point(0.721, 0.682), 0.465, 24, 10, false);
		auto g3 = GridGeom::cross_grids(&g1, &g2, 0.0, 0, false, false, 0, CrossGridCallback::silent());
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
		add_file_check(3784851216684557860, "g1.msh", "mesh with polyhedra cells");
		delete g3;
	}
}

int main(){
	test01();
	test02();
	test03();
	
	if (FAILED_CHECKS ==1){
		std::cout<<FAILED_CHECKS<<" test failed <<<<<<<<<<<<<<<<<<<"<<std::endl;
	} else if (FAILED_CHECKS > 1) {
		std::cout<<FAILED_CHECKS<<" tests failed <<<<<<<<<<<<<<<<<<<"<<std::endl;
	} else {
		std::cout<<"All tests passed"<<std::endl;
	}
	std::cout<<"DONE"<<std::endl;
}
