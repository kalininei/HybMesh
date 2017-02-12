#include "hmgrid3d.hpp"
#include <fstream>
#include "debug3d.hpp"
#include "hmtesting.hpp"
#include "hmtimer.hpp"
#include "buildgrid.hpp"
#include "contour.hpp"
#include "healgrid.hpp"
#include "assemble3d.hpp"

using namespace HMTesting;

void test01(){
	std::cout<<"1. Surface tree assembling, reverting, volumes"<<std::endl;
	{
		auto g1 = HM3D::Grid::Constructor::Cuboid({0, 0, 0}, 1, 1, 2, 2, 3, 2);
		std::swap(g1.vfaces[0], g1.vfaces[5]);
		auto s1 = HM3D::Surface::Assembler::GridSurface(g1);
		double v1, v2, v3, v4, v5, v6;
		v1 = HM3D::Surface::Volume(s1);
		{
			HM3D::Surface::R::Revert rr(s1);
			v2 = HM3D::Surface::Volume(s1);
			rr.reverse_direction();
			v3 = HM3D::Surface::Volume(s1);
			rr.reverse_direction();
			v4 = HM3D::Surface::Volume(s1);
			rr.reverse_direction();
			v5 = HM3D::Surface::Volume(s1);
			//now rr is deleted
		}
		v6 = HM3D::Surface::Volume(s1);
		add_check(!ISEQ(v1, v2) && ISEQ(v2, -2) && ISEQ(v3, 2) &&
		          ISEQ(v4, v2) && ISEQ(v5, v3) && ISEQ(v6, v1),
		          "cuboid surface temporal reverse procedure");
	}
	{
		auto gcyl2 = HM2D::Grid::Constructor::Circle(Point{0, 0}, 1, 64, 3, false);
		auto gcyl = HM3D::Grid::Constructor::SweepGrid2D(gcyl2, {0, 1, 2, 3});
		auto gtmp1 = HM2D::Grid::Constructor::Circle(Point(0, 0), 0.3, 64, 3, true);
		vector<int> inpcells;
		vector<int> badpoints;
		for (int i=0; i<gtmp1.vvert.size(); ++i){
			if (ISLOWER(gtmp1.vvert[i]->x, 0)) badpoints.push_back(i);
		}
		aa::enumerate_ids_pvec(gtmp1.vvert);
		for (int i=0; i<gtmp1.vcells.size(); ++i){
			auto c = gtmp1.vcells[i];
			bool good = true;
			auto op = HM2D::Contour::OrderedPoints(c->edges);
			for (int j=0; j<c->edges.size(); ++j){
				int pind = op[j]->id;
				if (std::find(badpoints.begin(), badpoints.end(), pind) !=
						badpoints.end()){
					good = false;
					break;
				}
			}
			if (good) inpcells.push_back(i);
		}
		HM2D::CellData cls2;
		for (auto i: inpcells) cls2.push_back(gtmp1.vcells[i]);
		auto ghsphere2 = HM2D::Grid::Constructor::InvokeGrid::Permanent(cls2);

		vector<double> degs = {0};
		for (int i=0; i<32; ++i) degs.push_back(degs.back() + 180./32.);
		auto ghsphere = HM3D::Grid::Constructor::RevolveGrid2D(ghsphere2, degs,
				Point(0, 0), Point(0, 1), false);
		for (auto p: ghsphere.vvert) p->z += 2.5;

		auto gcube1 = HM3D::Grid::Constructor::Cuboid({0, 0, 2.3}, 0.05, 0.05, 0.05, 2, 3, 4);
		auto gcube2 = HM3D::Grid::Constructor::Cuboid({20, 20, -2.5}, 6, 3, 1, 3, 1, 2);
		auto gcube3 = HM3D::Grid::Constructor::Cuboid({20, 20, -2.5}, 1, 1, 1, 3, 1, 2);

		HM3D::FaceData totalsurface;
		auto surf1 = HM3D::Surface::Assembler::GridSurface(ghsphere);
		auto surf2 = HM3D::Surface::Assembler::GridSurface(gcyl);
		auto surf3 = HM3D::Surface::Assembler::GridSurface(gcube1);
		auto surf4 = HM3D::Surface::Assembler::GridSurface(gcube2);
		auto surf5 = HM3D::Surface::Assembler::GridSurface(gcube3);
		totalsurface.insert(totalsurface.end(), surf1.begin(), surf1.end());
		totalsurface.insert(totalsurface.end(), surf2.begin(), surf2.end());
		totalsurface.insert(totalsurface.end(), surf3.begin(), surf3.end());
		totalsurface.insert(totalsurface.end(), surf4.begin(), surf4.end());
		totalsurface.insert(totalsurface.end(), surf5.begin(), surf5.end());
		std::random_shuffle(totalsurface.begin(), totalsurface.end());

		auto tree = HM3D::Surface::Tree::Assemble(totalsurface);
		double v1 = HM3D::Surface::Volume(totalsurface);
		auto* rr = new HM3D::Surface::R::RevertTree(tree);
		double v2 = HM3D::Surface::Volume(totalsurface);
		delete rr;
		double v3 = HM3D::Surface::Volume(totalsurface);
		add_check(fabs(v2 - 26.3534)<1e-4 && ISEQ(v1, v3), "complicated tree structure volume");
	}
}


int main(){
	test01();
	
	check_final_report();
	std::cout<<"DONE"<<std::endl;
}
