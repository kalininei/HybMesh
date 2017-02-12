#include "hmgrid3d.hpp"
#include <fstream>
#include "debug3d.hpp"
#include "hmtesting.hpp"
#include "hmtimer.hpp"
#include "buildgrid.hpp"
#include "unite_grids.hpp"
#include "healgrid.hpp"
#include "buildcont.hpp"
#include "trigrid.hpp"
#include "pebi.hpp"
#include "infogrid.hpp"
#include "assemble2d.hpp"
#include "assemble3d.hpp"
#include "debug2d.hpp"
#include "export2d_fluent.hpp"
#include "export2d_vtk.hpp"
using namespace HMTesting;

void old_numering(HM2D::GridData& g){
	aa::enumerate_ids_pvec(g.vvert);
	std::map<std::pair<int, int>, int> eds;
	for (int i=0; i<g.vedges.size(); ++i){
		auto ed = std::make_pair(g.vedges[i]->pfirst()->id, g.vedges[i]->plast()->id);
		if (ed.first > ed.second){
			std::swap(ed.first, ed.second);
			g.vedges[i]->reverse();
		}
		eds.emplace(ed, i);
	}
	HM2D::EdgeData newe(g.vedges.size());
	int ic=0;
	for (auto& it: eds){
		newe[ic++] = g.vedges[it.second];
	}
	std::swap(g.vedges, newe);
};

void test01(){
	std::cout<<"1. export cuboid to vtk"<<std::endl;
	auto g1 = HM3D::Grid::Constructor::Cuboid({0, 0, 0}, 1, 2, 5, 3, 3, 3);
	HM3D::Export::AllVTK.Silent(g1, "g1.vtk", "c1.vtk");
	add_check(g1.vvert.size() == 64 && g1.vcells.size() == 27 &&
			g1.vedges.size() == 144 && g1.vfaces.size() == 108,
			"cuboid primitives number"); 
	add_file_check(15732503264642486832U, "g1.vtk", "grid");
	add_file_check(12574868808136614456U, "c1.vtk", "boundary");
}

void test02(){
	std::cout<<"2. parallel sweep"<<std::endl;
	{
		auto g2d = HM2D::Grid::Constructor::RectGrid01(11, 7);
		auto g3d = HM3D::Grid::Constructor::SweepGrid2D(g2d, {0.3, 0.4, 0.8});
		HM3D::Export::BoundaryVTK(g3d, "c1.vtk");
		HM3D::Export::GridVTK(g3d, "g1.vtk");
		add_check(g3d.vvert.size() == 288 && g3d.vedges.size() == 708 &&
				g3d.vfaces.size() == 575 && g3d.vcells.size() == 154,
				"rectangular grid sweep");
	}
	{
		auto g2d = HM2D::Grid::Constructor::Ring(Point(0, 0), 4, 2, 12, 4);
		auto g3d = HM3D::Grid::Constructor::SweepGrid2D(g2d, {0, 0.1, 0.4, 0.7});
		HM3D::Export::BoundaryVTK(g3d, "c1.vtk");
		HM3D::Export::GridVTK(g3d, "g1.vtk");
		add_check(g3d.vvert.size() == 240 && g3d.vcells.size() == 144,
				"ring grid sweep");
	}
	{
		auto g2d = HM2D::Grid::Constructor::Circle(Point(1, 0), 4, 24, 10, true);
		auto g3d = HM3D::Grid::Constructor::SweepGrid2D(g2d, {0, 1, 2, 3});
		HM3D::Export::BoundaryVTK(g3d, "c1.vtk");
		HM3D::Export::GridVTK(g3d, "g1.vtk");
		add_check(g3d.vvert.size() == 964 && g3d.vcells.size() == 720,
				"circle grid sweep");
	}
}

void test03(){
	std::cout<<"3. Fluent export"<<std::endl;
	{
		auto g1 = HM3D::Grid::Constructor::Cuboid(HM3D::Vertex(0, 0, 0), 1, 1, 1, 2, 2, 1);
		HM3D::Export::GridMSH.Silent(g1, "g1.msh");
		add_file_check(15306802383608290446U, "g1.msh", "simple cuboid");
	}
	{
		auto g2d = HM2D::Grid::Constructor::RectGrid01(6, 3);
		old_numering(g2d);
		auto g3d = HM3D::Grid::Constructor::SweepGrid2D(g2d, {0, 0.1, 0.2, 0.5},
				[](int i){ return 1; },
				[](int i){ return 2; },
				[](int i){ return i+3; });
		HM3D::Export::GridMSH(g3d, "g2.msh",
				[](int i)->std::string{
					switch (i){
						case 1: return std::string("bottom");
						case 2: return std::string("top");
						default: return std::string("side") + std::to_string(i);
					};
				});
		add_file_check(7968177351678915047U, "g2.msh", "cuboid from sweep with custom boundaries");
	}
	{
		auto g1 = HM2D::Grid::Constructor::Circle(Point(0, 0), 1, 4, 2, true);
		HM2D::CellData cd {g1.vcells[0], g1.vcells[4]};
		auto g2 = HM2D::Grid::Constructor::InvokeGrid::Permanent(cd);
		old_numering(g2);
		auto g3d = HM3D::Grid::Constructor::SweepGrid2D(g2, {0, 0.1});
		HM3D::Export::GridMSH(g3d, "g3.msh");
		add_file_check(3132562703476878584U, "g3.msh", "mixed hex/wedge cells");
	}
	{
		auto g1 = HM2D::Grid::Constructor::Circle(Point(0, 0), 1, 5, 2, false);
		HM2D::CellData cd {g1.vcells[5]};
		auto g2 = HM2D::Grid::Constructor::InvokeGrid::Permanent(cd);
		old_numering(g2);
		auto g3d = HM3D::Grid::Constructor::SweepGrid2D(g2, {0, 0.1});
		HM3D::Export::GridMSH(g3d, "g4.msh");
		add_file_check(1811066807055341946U, "g4.msh", "single pentagon prism cell");
	}
	{
		auto g1 = HM2D::Grid::Constructor::RectGrid01(20, 30);
		old_numering(g1);
		auto g2 = HM2D::Grid::Constructor::Circle(Point(0.721, 0.682), 0.465, 24, 10, false);
		old_numering(g2);
		auto g3 = HM2D::Grid::Algos::UniteGrids(g1, g2, HM2D::Grid::Algos::OptUnite());
		old_numering(g3);
		auto g3d = HM3D::Grid::Constructor::SweepGrid2D(g3, {0, 0.1, 0.2, 0.3, 0.5},
				[](int){return 1;}, [](int){return 2;},
				[&g3](int e)->int{
					Point pc = g3.vedges[e]->center();
					return (pc.x<=1+1e-12 && pc.y<=1+1e-12) ? 3 : 4;
				});
		HM3D::Export::GridMSH(g3d, "g5.msh",
				[](int i)->std::string{
					switch (i){
						case 1: return "bottom";
						case 2: return "top";
						case 3: return "square";
						case 4: return "circle";
						default: return "unknown";
					}
				});
		add_file_check(4247503388199499266U, "g5.msh", "mesh with polyhedra cells");
	}
}

void test04(){
	std::cout<<"4. Fluent export with periodic surfaces"<<std::endl;
	{
		auto g2d = HM2D::Grid::Constructor::RectGrid({0.0, 0.1, 1.0}, {0.0, 0.3, 1.0});
		old_numering(g2d);
		auto g3d = HM3D::Grid::Constructor::SweepGrid2D(g2d, {0.0, 0.2, 1.0});
		HM3D::Export::PeriodicData pd;
		pd.add_condition(1, 2, HM3D::Vertex(0, 0, 0), HM3D::Vertex(0, 0, 1), true);
		HM3D::Export::GridMSH(g3d, "_o1.msh", pd); 
		add_file_check(1901761016274060527U, "_o1.msh", "simple 2x2x2");

		pd.data[0].reversed = false;
		HM3D::Export::GridMSH(g3d, "_o2.msh", pd); 
		add_file_check(17909037251898648897U, "_o2.msh", "2x2x2 without reverse");

		pd.data[0].reversed = true;
		pd.data[0].v = HM3D::Vertex(0.1, 0, 0);
		bool res3 = false;
		try{
			HM3D::Export::GridMSH(g3d, "g3.msh", pd); 
		} catch(std::runtime_error& e){
			res3 = true;
		}
		add_check(res3, "fail at invalid point match");
	}
	{
		auto g2d1 = HM2D::Grid::Constructor::RectGrid(Point(0, 0), Point(10, 1), 100, 10);
		auto g2d2 = HM2D::Grid::Constructor::Ring(Point(3, 0.5), 0.3, 0.1, 20, 4);
		HM2D::Grid::Algos::OptUnite opt(0.1);
		opt.empty_holes = true;
		auto g2d = HM2D::Grid::Algos::UniteGrids(g2d1, g2d2, opt);
		vector<double> zvec;
		for (int i=0; i<100; i+=10)  zvec.push_back(3 + (double)i/99);
		auto g3d = HM3D::Grid::Constructor::SweepGrid2D(g2d, zvec);

		HM3D::Ser::Grid s3d(g3d);
		s3d.set_btype([](HM3D::Vertex v, int bt){
					if (bt == 3){
						if (v.x<=geps) return 3;
						if (v.x>=10-geps) return 4;
						return 5;
					}
					return bt;
				});

		HM3D::Export::PeriodicData pd;

		pd.add_condition(1, 2, HM3D::Vertex(0, 0, 3), HM3D::Vertex(0, 0, 4), true);
		pd.add_condition(3, 4, HM3D::Vertex(0, 0, 3), HM3D::Vertex(10, 0, 3), true);
		HM3D::Export::GridMSH(g3d, "g2.msh", pd);
		add_file_check(18279916701421103105U, "g2.msh", "multiple periodic");
	}
};

void test05(){
	std::cout<<"5. Tecplot export"<<std::endl;
	{
		auto g2d = HM2D::Grid::Constructor::RectGrid01(1, 1);
		old_numering(g2d);
		auto g3d = HM3D::Grid::Constructor::SweepGrid2D(g2d, {0.0, 0.5});
		HM3D::Export::GridTecplot.Silent(g3d, "g1.dat");
		add_file_check(1831833575709478659U, "g1.dat", "single cell grid");
	}
	{
		auto g2d = HM2D::Grid::Constructor::Circle(Point(0, 0), 10, 30, 10, false);
		old_numering(g2d);
		auto g3d = HM3D::Grid::Constructor::SweepGrid2D(g2d, {1.0, 1.2, 1.4, 1.6, 1.7, 1.8, 1.9, 2.0});
		HM3D::Export::GridTecplot.Silent(g3d, "g1.dat");
		add_file_check(17626851046985520587U, "g1.dat", "polyhedral grid");
		HM3D::Export::BoundaryTecplot.Silent(HM3D::Ser::Grid(g3d), "g1.dat");
		add_file_check(8291026423155100327U, "g1.dat", "polyhedral boundary");
	}
}

void test06(){
	using HM3D::Grid::Constructor::RevolveGrid2D;
	std::cout<<"6. Solid of revolution"<<std::endl;
	auto g2d = HM2D::Grid::Constructor::RectGrid(Point(1,0), Point(2,1), 1, 1);
	old_numering(g2d);
	auto bc0 = [](int){return 0;};
	{
		for (auto e: g2d.vedges) e->boundary_type = 0;
		HM3D::GridData g3d = RevolveGrid2D(g2d, {0, 90}, Point(0, 0), Point(0, 1),
				true, 0, 0);
		HM3D::Export::GridTecplot.Silent(g3d, "g1.dat");
		add_file_check(16088294825526263046U, "g1.dat", "single cell, distant, incomplete");
	}
	{
		for (auto e: g2d.vedges) e->boundary_type = 0;
		auto g3d = RevolveGrid2D(g2d, {0, 90, 180, 270, 360}, Point(0, 0), Point(0, 1),
				true, 0, 0);
		HM3D::Export::GridTecplot.Silent(g3d, "g1.dat");
		add_file_check(8732440237994901672U, "g1.dat", "single cell, distant, complete");
	}
	{
		for (int i=0; i<g2d.vedges.size(); ++i) g2d.vedges[i]->boundary_type = i;
		auto g3d = RevolveGrid2D(g2d, {0, 90, 100}, Point(0, 0), Point(0, 1),
				true, 10, 20);
		HM3D::Export::GridTecplot.Silent(g3d, "g1.dat");
		add_file_check(3859847262675033285U, "g1.dat", "single cell, distant, incomplete, with bc");
	}
	{
		auto h2d = HM2D::Grid::Constructor::RectGrid(Point(0, 0), Point(2, 1), 2, 1);
		old_numering(h2d);
		for (auto& e: h2d.vedges) e->boundary_type = 1;
		auto g3d = RevolveGrid2D(h2d, {0, 90}, Point(0, 0), Point(0, 1), true);
		HM3D::Export::GridTecplot.Silent(g3d, "g1.dat");
		add_file_check(8233442907809870919U, "g1.dat", "with contact, incomplete");
	}
	{
		auto h2d = HM2D::Grid::Constructor::RectGrid(Point(0, 0), Point(2, 1), 4, 3);
		old_numering(h2d);
		for (auto& e: h2d.vedges) e->boundary_type = 1;
		auto g3d = RevolveGrid2D(h2d, {0, 90, 110, 180, 250, 330, 360}, Point(0, 0), Point(0, 1), true);
		HM3D::Export::GridTecplot.Silent(g3d, "g1.dat");
		add_file_check(5490115627065179709U, "g1.dat", "with contact, complete");
	}
	{
		auto g1 = HM2D::Grid::Constructor::RectGrid(Point(0, 0), Point(10, 10), 10, 10);
		old_numering(g1);
		auto g2 = HM2D::Grid::Constructor::RectGrid(Point(0, 5), Point(10, 6), 5, 1);
		old_numering(g2);
		auto g3 = HM2D::Grid::Algos::UniteGrids(g1, g2, HM2D::Grid::Algos::OptUnite());
		old_numering(g3);
		for (auto& e: g3.vedges) e->boundary_type = 1;
		auto g3d = RevolveGrid2D(g3, {0, 10, 20, 30}, Point(0, 0), Point(0, 1), true);
		HM3D::Export::GridTecplot.Silent(g3d, "g1.dat");
		add_file_check(12980710001405184230U, "g1.dat", "hanging nodes near axis to tecplot");
		HM3D::Export::GridMSH.Silent(g3d, "g1.msh");
		add_file_check(8061023987823183823U, "g1.msh", "hanging nodes near axis to fluent");
	}
}

void test07(){
	using HM3D::Grid::Constructor::RevolveGrid2D;
	std::cout<<"7. Solid of revolution, merging centeral cells"<<std::endl;
	{
		auto g2d = HM2D::Grid::Constructor::RectGrid(Point(1,0), Point(2,1), 1, 1);
		old_numering(g2d);
		for (auto e: HM2D::ECol::Assembler::GridBoundary(g2d)) e->boundary_type = 1;
		auto g3d = RevolveGrid2D(g2d, {0, 45, 90}, Point(1, 0), Point(1, 1), false);
		HM3D::Export::GridTecplot.Silent(g3d, "g1.dat");
		add_file_check(13398422286724743124U, "g1.dat", "single cell, without center trian, incomplete");
	}
	{
		auto g2d = HM2D::Grid::Constructor::RectGrid(Point(1,0), Point(2,1), 1, 1);
		for (auto e: HM2D::ECol::Assembler::GridBoundary(g2d)) e->boundary_type = 1;
		old_numering(g2d);
		auto g3d = RevolveGrid2D(g2d, {20, 45, 90, 160, 270, 300, 380}, Point(1, 0), Point(1, 1), false);
		HM3D::Export::GridTecplot.Silent(g3d, "g1.dat");
		add_file_check(6994418583934313116U, "g1.dat", "single cell, without center trian, complete");
	}
	{
		auto h2d = HM2D::Grid::Constructor::RectGrid(Point(0, 0), Point(2, 1), 2, 1);
		for (auto e: HM2D::ECol::Assembler::GridBoundary(h2d)) e->boundary_type = 1;
		old_numering(h2d);
		auto g3d = RevolveGrid2D(h2d, {0, 10, 20, 30, 40, 50}, Point(0, 0), Point(0, 1), false);
		HM3D::Export::GridTecplot.Silent(g3d, "g1.dat");
		add_file_check(11881236001573517783U, "g1.dat", "multiple cells, with trian, complete");
	}
	{
		double p[] = {0, 0, 1, 0, 0, 1};
		int c[] = {0, 1, 2};
		auto g1 = HM2D::Grid::Constructor::FromRaw(3, 1, p, c, 3);
		old_numering(g1);
		for (auto e: HM2D::ECol::Assembler::GridBoundary(g1)) e->boundary_type = 1;
		auto g2 = RevolveGrid2D(g1, {0, 45, 90}, Point(0,0), Point(0,1), false);
		HM3D::Export::GridTecplot.Silent(g2, "g1.dat");
		add_file_check(10167032458429566145U, "g1.dat", "no tri with single axis triangle");
	}
	{
		double p[] = {1, 0, 1, 1, 0, 1};
		int c[] = {0, 1, 2};
		auto g1 = HM2D::Grid::Constructor::FromRaw(3, 1, p, c, 3);
		old_numering(g1);
		for (auto e: HM2D::ECol::Assembler::GridBoundary(g1)) e->boundary_type = 1;
		auto g2 = RevolveGrid2D(g1, {0, 45, 90}, Point(0,0), Point(0,1), false);
		HM3D::Export::GridTecplot.Silent(g2, "g1.dat");
		add_file_check(11550191908304285294U, "g1.dat", "no tri, single off axis triangle");
	}
	{
		double p[] = {0, 0,  1, 0,  0, 1,
		              1, 0,  1, 1,  0, 1,
		              0, 0,  0, -2,  1, -2,  1, 0};
		int c[] = {3, 0, 1, 2, 
		           3, 3, 4, 5,
		           4, 6, 7, 8, 9};
		auto g1 = HM2D::Grid::Constructor::FromRaw(10, 3, p, c, -1);
		old_numering(g1);
		HM2D::Grid::Algos::Heal(g1);
		old_numering(g1);
		for (auto e: HM2D::ECol::Assembler::GridBoundary(g1)) e->boundary_type = 1;
		auto g2 = RevolveGrid2D(g1, {0, 45, 90}, Point(0,0), Point(0,1), false);
		HM3D::Export::GridTecplot.Silent(g2, "g1.dat");
		add_file_check(12664340621499564857U, "g1.dat", "no tri, complex connections, incomplete");
		auto g3 = RevolveGrid2D(g1, {0, 90, 180, 270, 360}, Point(0,0), Point(0,1), false);
		HM3D::Export::GridTecplot.Silent(g3, "g1.dat");
		add_file_check(2848618625331037303U, "g1.dat", "no tri, complex connections, complete");
	}
}

void test08(){
	std::cout<<"8. export cuboid to gmsh"<<std::endl;
	auto g1 = HM3D::Grid::Constructor::Cuboid({0, 0, 0}, 1, 2, 5, 3, 3, 3);
	auto bfun = [](int i)->std::string{
		if (i==1) return "bottom";
		if (i==2) return "top";
		if (i==3) return "left";
		if (i==4) return "right";
		if (i==5) return "front";
		return "back";
	};
	HM3D::Export::GridGMSH(HM3D::Ser::Grid(g1), "g1.msh", bfun);
	add_check(g1.vvert.size() == 64 && g1.vcells.size() == 27 &&
			g1.vedges.size() == 144 && g1.vfaces.size() == 108,
			"cuboid primitives number"); 
	add_file_check(4596785021162173517U, "g1.msh", "3d gmsh export");
}

void test09(){
	std::cout<<"9. 3d domain unstructured meshing"<<std::endl;
	{
		auto g1 = HM3D::Grid::Constructor::Cuboid({0, 0, 0}, 1, 1, 1, 5, 5, 5);
		auto s1 = HM3D::Surface::Assembler::GridSurface(g1);
		auto g2 = HM3D::Mesher::UnstructuredTetrahedral(s1); 
		add_check(ISEQ(HM3D::SumVolumes(g2.vcells), 1), "grid in cubic domain");
	}
	{
		auto g1 = HM3D::Grid::Constructor::Cuboid({1, 1, 1}, 2, 3, 1, 7, 8, 4);
		auto g2 = HM3D::Grid::Constructor::Cuboid({10, 10, 10}, 5, 5, 5, 10, 10, 10);
		auto gcyl2 = HM2D::Grid::Constructor::Circle(Point{1, 1}, 5, 64, 10, true);
		vector<double> zsweep;
		for (int i=0; i<11; ++i) zsweep.push_back( (double)i / 10 * 4);
		auto gcyl = HM3D::Grid::Constructor::SweepGrid2D(gcyl2, zsweep);

		HM3D::FaceData srf;
		//for (auto f: gcyl.vfaces) if (f->is_boundary()) srf.faces.push_back(f);
		//for (auto f: g1.vfaces) if (f->is_boundary()) srf.faces.push_back(f);
		for (auto f: g2.vfaces) if (f->is_boundary()) srf.push_back(f);

		auto res = HM3D::Mesher::UnstructuredTetrahedral.WVerbTimer(srf);
		//auto res = HM3D::Mesher::UnstructuredTetrahedral(srf);
		HM3D::Export::SurfaceVTK(srf, "srf.vtk");
		HM3D::Export::GridVTK(res, "res.vtk");

		double v1 = HM3D::SumVolumes(res.vcells);
		auto tree = HM3D::Surface::Tree::Assemble(srf);
		auto* rr = new HM3D::Surface::R::RevertTree(tree);
		double v2 = HM3D::Surface::Volume(srf);
		delete rr;
		add_check(ISEQ(v1, v2), "multiply connected domain");
	}
	{
		vector<Point> pts = {Point(0, 0), Point(1, 0), Point(1.2, 1), Point(-0.1, 0.9)};
		HM2D::EdgeData cont = HM2D::Contour::Constructor::FromPoints(pts, true);
		auto tree = HM2D::Mesher::PrepareSource(cont, 0.1);
		auto g1 = HM2D::Mesher::UnstructuredTriangle(tree);
		auto g2 = HM2D::Grid::Constructor::TriToPebi(g1);
		vector<double> zsweep;
		for (int i=0; i<11; ++i) zsweep.push_back( (double)i / 10 );
		auto g3 = HM3D::Grid::Constructor::SweepGrid2D(g2, zsweep);

		HM3D::FaceData srf;
		for (auto f: g3.vfaces) if (f->is_boundary()) srf.push_back(f);
		auto res = HM3D::Mesher::UnstructuredTetrahedral.WVerbTimer(srf);

		double v1 = HM2D::Grid::Area(g1) * 1.;
		add_check(ISEQ(v1, HM3D::SumVolumes(res.vcells)), "pebi base");
	}
	{
		vector<Point> pts = {Point(0, 0), Point(1.2, 0.3), Point(1.1, 1), Point(-0.3, 0.9)};
		HM2D::EdgeData cont = HM2D::Contour::Constructor::FromPoints(pts, true);
		auto tree = HM2D::Mesher::PrepareSource(cont, 0.1);
		auto g1 = HM2D::Mesher::UnstructuredTriangle(tree);
		vector<double> zsweep;
		for (int i=0; i<11; ++i) zsweep.push_back( (double)i / 10 );
		auto g3 = HM3D::Grid::Constructor::SweepGrid2D(g1, zsweep);

		HM3D::FaceData srf;
		for (auto f: g3.vfaces) if (f->is_boundary()) srf.push_back(f);
		auto res = HM3D::Mesher::UnstructuredTetrahedral.WVerbTimer(srf);

		double v1 = HM2D::Grid::Area(g1) * 1.;
		add_check(ISEQ(v1, HM3D::SumVolumes(res.vcells)), "triangulated base");
		HM3D::Export::GridVTK(res, "g1.vtk");
	}
	{
		vector<Point> pts = {Point(0, 0), Point(1.2, 0.3), Point(1.1, 1)};
		HM2D::EdgeData cont = HM2D::Contour::Constructor::FromPoints(pts, true);
		auto tree = HM2D::Mesher::PrepareSource(cont, 0.1);
		auto g1 = HM2D::Mesher::UnstructuredTriangle(tree);
		vector<double> zsweep;
		for (int i=0; i<11; ++i) zsweep.push_back( (double)i / 10 );
		auto g3 = HM3D::Grid::Constructor::SweepGrid2D(g1, zsweep);

		HM3D::FaceData srf;
		for (auto f: g3.vfaces) if (f->is_boundary()) srf.push_back(f);
		auto res = HM3D::Mesher::UnstructuredTetrahedral.WVerbTimer(srf);

		double v1 = HM2D::Grid::Area(g1) * 1.;
		add_check(ISEQ(v1, HM3D::SumVolumes(res.vcells)), "acute angle in the area");
		HM3D::Export::GridVTK(res, "g1.vtk");
	}
	{
		HM2D::GridData g1 = HM2D::Grid::Constructor::RegularHexagonal(Point(0, 0), Point(10, 10), 1);
		vector<double> zsweep {0};
		//FIXME hz=0.2 fails due to pyramid intersections
		double hz = 0.25;
		while (zsweep.back()<5) zsweep.push_back(zsweep.back()+hz);
		auto g2 = HM3D::Grid::Constructor::SweepGrid2D(g1, zsweep);
		HM3D::FaceData srf2;
		for (auto f: g2.vfaces) if (f->is_boundary()) srf2.push_back(f);

		auto res2 = HM3D::Mesher::UnstructuredTetrahedral.WVerbTimer(srf2);
		add_check(ISEQ(HM2D::Grid::Area(g1)*zsweep.back(),
		               HM3D::SumVolumes(res2.vcells)),
		          "non-smooth boundaries");
	}
}
void test10(){
	std::cout<<"10. Meshing of volumes with bottlenecks"<<std::endl;
	HM2D::GridData g1 = HM2D::Grid::Constructor::RegularHexagonal(Point(0, 0), Point(10, 10), 1);
	HM2D::CellData cd {g1.vcells[0]};
	HM2D::GridData g2 = HM2D::Grid::Constructor::InvokeGrid::Permanent(cd);
	{
		vector<double> zsweep;
		double hz = 1;
		for (int i=0; i<2; ++i) zsweep.push_back(i*hz);
		auto g3 = HM3D::Grid::Constructor::SweepGrid2D(g2, zsweep);
		HM3D::FaceData srf3;
		for (auto f: g3.vfaces) if (f->is_boundary()) srf3.push_back(f);

		auto res3 = HM3D::Mesher::UnstructuredTetrahedral(srf3);
		add_check(ISEQ(HM2D::Grid::Area(g2)*hz, HM3D::SumVolumes(res3.vcells)),
				"1 layer: whole area filled with pyramids");
	}
	{
		vector<double> zsweep;
		double hz = 1;
		for (int i=0; i<3; ++i) zsweep.push_back(i*hz);
		auto g3 = HM3D::Grid::Constructor::SweepGrid2D(g2, zsweep);
		HM3D::FaceData srf3;
		for (auto f: g3.vfaces) if (f->is_boundary()) srf3.push_back(f);

		auto res3 = HM3D::Mesher::UnstructuredTetrahedral(srf3);
		add_check(ISEQ(HM2D::Grid::Area(g2)*2*hz, HM3D::SumVolumes(res3.vcells)),
				"2 layers");
	}
	{
		//FIXME: this test doesn't work because built
		//pyramids change nesting structure
		/*
		vector<double> zsweep;
		double hz = 1;
		for (int i=0; i<4; ++i) zsweep.push_back(i*hz);
		auto g3 = HM3D::Grid::Constructor::SweepGrid2D(g2, zsweep);
		HM3D::Surface srf3;
		for (auto f: g3.vfaces) if (f->is_boundary()) srf3.faces.push_back(f);

		auto res3 = HM3D::Mesher::UnstructuredTetrahedral.WVerbTimer(srf3);
		add_check(ISEQ(HM2D::Grid::Info::Area(g2)*3*hz, HM3D::Cell::SumVolumes(res3.vcells)),
				"3 layers");
		*/
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
	test08();
	test09();
	test10();
	
	check_final_report();
	std::cout<<"DONE"<<std::endl;
}
