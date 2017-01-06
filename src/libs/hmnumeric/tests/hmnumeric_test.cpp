#include "hmfem.hpp"
#include "femassembly.hpp"
#include "hmtesting.hpp"
#include "vtk_export_grid2d.hpp"
#include "cont_assembler.hpp"
#include "constructor.hpp"
#include "vtk_export2d.hpp"
#include "cont_partition.hpp"
using HMTesting::add_check;

void test01(){
	std::cout<<"01. fem in square"<<std::endl;
	Point p1{0.0, 0.0}, p2{2.0, 0.0}, p3{2.0, 1.0}, p4{0.0, 1.0};

	//-> Dfdn(bottom), Dfdx(area), DfDy(area)
	auto dotest = [&](shared_ptr<HMFem::Grid43> g)->std::array<double, 3>{
		auto cont = GGeom::Info::Contour1(*g);
		auto _av = HM2D::AllVertices(cont);
		auto cp1 = _av[std::get<0>(HM2D::FindClosestNode(_av, p1))].get();
		auto cp2 = _av[std::get<0>(HM2D::FindClosestNode(_av, p2))].get();
		auto cp3 = _av[std::get<0>(HM2D::FindClosestNode(_av, p3))].get();
		auto cp4 = _av[std::get<0>(HM2D::FindClosestNode(_av, p4))].get();

		//GridPoint* gp1 = static_cast<GridPoint*>(HMCont2D::ECollection::FindClosestNode(cont, p1));
		//GridPoint* gp2 = static_cast<GridPoint*>(HMCont2D::ECollection::FindClosestNode(cont, p2));
		//GridPoint* gp3 = static_cast<GridPoint*>(HMCont2D::ECollection::FindClosestNode(cont, p3));
		//GridPoint* gp4 = static_cast<GridPoint*>(HMCont2D::ECollection::FindClosestNode(cont, p4));
		auto contbot = HM2D::Contour::Assembler::ShrinkContour(cont, cp1, cp2);
		auto conttop = HM2D::Contour::Assembler::ShrinkContour(cont, cp3, cp4);

		auto p = HMFem::LaplasProblem(g);
		p.SetDirichlet(contbot, [](const GridPoint* p){ return 1.0; });
		p.SetDirichlet(conttop, [](const GridPoint* p){ return 6.0; });


		auto v = vector<double>(g->n_points(), 0.0);
		p.Solve(v);

		double I1 = p.IntegralDfDn(contbot, v);

		auto mass = HMFem::Assemble::LumpMass(*g);

		vector<double> dx(g->n_points());
		HMFem::Assemble::DDx(*g)->MultVec(v, dx);
		for (int i=0; i<dx.size(); ++i) dx[i]/=mass[i];
		double I2 = std::inner_product(dx.begin(), dx.end(), mass.begin(), 0.0);

		vector<double> dy(g->n_points());
		HMFem::Assemble::DDy(*g)->MultVec(v, dy);
		for (int i=0; i<dy.size(); ++i) dy[i]/=mass[i];
		double I3 = std::inner_product(dy.begin(), dy.end(), mass.begin(), 0.0);

		return {I1, I2, I3};
	};

	//1) structured grid
	auto g4 = GGeom::Constructor::RectGrid01(21, 21);
	GGeom::Modify::PointModify(g4, [](GridPoint* p){ p->x*=2; });
	auto g1 = HMFem::Grid43::Build3(&g4);
	auto I1 = dotest(g1);
	add_check(fabs(I1[0]+10)<1e-6, "structured triangle grid, DfDn");
	add_check(fabs(I1[1])<1e-6, "structured triangle grid, DfDx");
	add_check(fabs(I1[2]-10)<1e-6, "structured triangle grid, DfDy");

	//2) unstructured grid
	double h = 0.05;
	auto g2 = HMFem::Grid43::Build3( {p1, p2, p3, p4}, h);
	auto I2 = dotest(g2);
	add_check(fabs(I2[0]+10)<1e-6, "unstructured triangle grid, DfDn");
	add_check(fabs(I2[1])<1e-6, "unstructured triangle grid, DfDx");
	add_check(fabs(I2[2]-10)<1e-6, "unstructured triangle grid, DfDy");
}

void test02(){
	std::cout<<"02. Auxiliary triangulation"<<std::endl;
	{
		//Rectangle checks
		auto rcont = HM2D::Contour::Constructor::FromPoints({0,0, 10,0, 10,10, 0,10}, true);
		auto ans1 = HMFem::AuxGrid3(rcont, 200, 200000);
		add_check(abs(ans1.n_points()-200)<0.1*200, "coarse grid in rectangle");

		auto ans3 = HMFem::AuxGrid3(rcont, 10000, 200000);
		add_check(abs(ans3.n_points()-10000)<0.1*10000, "fine grid in rectangle");
	}
	{
		//ring checks
		auto rgrid = GGeom::Constructor::Ring(Point(0, 0), 1, 0.7, 16, 1);
		auto rcont = GGeom::Info::Contour(rgrid);
		GridGeom ans1 = HMFem::AuxGrid3(rcont, 5000, 1000000);
		add_check(abs(ans1.n_points()-5000)<0.1*5000, "grid in the ring area");
	}
	{
		//grid with coarsening needed
		auto rcont = HM2D::Contour::Constructor::Circle(400, 1.0, Point(0, 0));
		GridGeom ans1 = HMFem::AuxGrid3(rcont, 1000, 1000000);
		add_check(abs(ans1.n_points()-1000)<0.5*1000, "grid in fine circle contour");
	}
	{
		auto rcont2 = HM2D::Contour::Constructor::FromPoints({
			0.0,0.0, 3.0,0.0,   3.0,1.0,   2.5,1.0, 2.47,1.04, 2.43,1,
			1.6,1.0, 1.57,1.04, 1.55,0.99, 0.0,1.0}, true);
		std::map<double, double> w;
		w[0] = 0.02; w[0.5] = 0.4;
		HM2D::Export::ContourVTK(rcont2, "c1.vtk");
		auto rcont = HM2D::Contour::Algos::WeightedPartition(w, rcont2,
				HM2D::Contour::Algos::PartitionTp::KEEP_ALL);
		HM2D::Export::ContourVTK(rcont, "c1.vtk");
		GridGeom ans1 = HMFem::AuxGrid3(rcont, 500, 1000000);
		GGeom::Export::GridVTK(ans1, "g1.vtk");
		add_check(abs(ans1.n_points()-500)<500, "grid with fine regions in its source tree");
	}
}

void test03(){
	std::cout<<"03. Auxiliary triangulation with contour constrains"<<std::endl;
	{
		auto sqrcont = HM2D::Contour::Constructor::Circle(32, 1.0, Point(0, 0));
		auto constr = HM2D::Contour::Constructor::FromPoints({
				-0.7,0, -0.69,0, -0.68,0, -0.67,0, -0.3,0, 0.1,0.4, 0.6,-0.3});
		HM2D::Contour::Tree tree;
		tree.AddContour(sqrcont);
		GridGeom ans1 = HMFem::AuxGrid3(tree, vector<HM2D::EdgeData> {constr}, 1000, 1000000);
		add_check([&](){
			auto fpoint = [&](double x, double y){
				Point pp(x, y);
				for (int i=0; i<ans1.n_points(); ++i)
					if (*ans1.get_point(i) == pp) return true;
				return false;
			};
			return (fpoint(-0.7, 0) && fpoint(-0.69, 0) && fpoint(-0.68,0) && fpoint(-0.67,0)
				&& fpoint(-0.3, 0) && ans1.n_points() < 1200);
		}(), "non-touching constraint");
	}
	{
		auto sqrcont = HM2D::Contour::Constructor::FromPoints({0,0, 1,0, 1,1, 0,1}, true);
		std::map<double, double> w;
		w[0.133] = 0.01; w[0.5] = 0.1;
		auto sqrcont2 = HM2D::Contour::Algos::WeightedPartition(w, sqrcont);
		auto lineconst = HM2D::Contour::Constructor::FromPoints({0.5,0, 0.2,0.5});
		HM2D::Contour::Tree tree;
		tree.AddContour(sqrcont2);
		GridGeom ans1 = HMFem::AuxGrid3(tree,
				vector<HM2D::EdgeData> {lineconst}, 500, 100000);
		add_check(ans1.n_points() < 700 && ans1.n_cells() < 1200 &&
			[&](){
				Point p1(0.5, 0), p2(0.2, 0.5);
				double ksieta[2];
				for (auto e: ans1.get_edges()){
					Point p3 = *ans1.get_point(e.p1);
					Point p4 = *ans1.get_point(e.p2);
					if (SectCross(p1, p2, p3, p4, ksieta)){
						if (ksieta[1]>geps && ksieta[1]<1-geps) return false;
					}
				}
				return true;
			}()
			, "touching coarse constraint");
	}
	{
		auto sqrcont = HM2D::Contour::Constructor::FromPoints({0,0, 1,0, 1,1, 0,1}, true);
		std::map<double, double> w;
		w[0.133] = 0.01; w[0.5] = 0.2;
		auto sqrcont2 = HM2D::Contour::Algos::WeightedPartition(w, sqrcont);
		auto lineconst = HM2D::Contour::Constructor::FromPoints({0.5,0, 0.2,0.5, 1,0.2});
		auto lineconst2 = HM2D::Contour::Algos::Partition(0.01, lineconst, HM2D::Contour::Algos::PartitionTp::KEEP_ALL);
		HM2D::Contour::Tree tree;
		tree.AddContour(sqrcont2);
		GridGeom ans1 = HMFem::AuxGrid3(tree,
				vector<HM2D::EdgeData> {lineconst2}, 500, 100000);
		add_check(ans1.n_points() < 800, "touching fine constraint");
	}
};


int main(){
	test01();
	test02();
	test03();

	HMTesting::check_final_report();
	std::cout<<"DONE"<<std::endl;
}
