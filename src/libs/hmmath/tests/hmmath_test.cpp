#include "fileproc.h"
#include "hybmesh_contours2d.hpp"
#include "scpack_port.hpp"
#include "dscpack_port.hpp"
#include "hmfem.hpp"
#include "fileproc.h"
#include "confrect_fem.hpp"
#include "femassembly.hpp"

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

void test01(){
	using namespace HMMath::Conformal::Impl::SCPack;
	std::cout<<"SCPACK conformal mapping"<<std::endl;
	auto left = HMCont2D::Constructor::ContourFromPoints({
			Point {0.4, 0.2},
			Point {0.47662610, 0.22298783},
		});
	auto bottom = HMCont2D::Constructor::ContourFromPoints({
			Point {0.4, 0.2},
			Point {0.46, 0.0},
		});
	auto right = HMCont2D::Constructor::ContourFromPoints({
			Point {0.46, 0.0},
			Point {0.53990000, -0.00000002},
			});
	auto top = HMCont2D::Constructor::ContourFromPoints({
			Point {0.47662610, 0.22298783},
			Point {0.53662610, 0.02298783},
			Point {0.53909957, 0.01196903},
			Point {0.53989935, 0.00401162},
			Point {0.53990000, -0.00000002},
		});
	auto trans = ToRect::Build(left, right, bottom, top);
	add_check(fabs(trans->module() - 2.73804)<1e-3, "Module");
	auto rect = trans->RectContour();
	double w = 0.678678;
	vector<Point> rectpnt { Point(w*trans->module(), 0), Point(w*trans->module(), 1) };
	vector<Point> polypnt = trans->MapToPolygon(rectpnt);
	add_check(
		fabs(polypnt[0].x - 0.442541 )<1e-6 &&
		fabs(polypnt[0].y - 0.0581967)<1e-6 &&
		fabs(polypnt[1].x - 0.519522 )<1e-6 &&
		fabs(polypnt[1].y - 0.0800001)<1e-6, 
		"Map to Polygon");
}

void test02(){
	using namespace HMMath::Conformal::Impl;
	std::cout<<"Rectangle approximation of conformal mapping"<<std::endl;
	auto top = HMCont2D::Constructor::ContourFromPoints({
		Point{0.43701993, 0.07659573},
		Point{0.44410482, 0.07840499},
		Point{0.45201250, 0.07960025},
		Point{0.46000000, 0.08000000},
		Point{1.00000000, 0.08000000}});

	auto right = HMCont2D::Constructor::ContourFromPoints({
		Point(1,0), Point(1,0.08)});

	auto bottom = HMCont2D::Constructor::ContourFromPoints({
		Point(0.46, 0.0), Point(1.0, 0.0)});

	auto left = HMCont2D::Constructor::ContourFromPoints({
		Point(0.46, 0.0), Point(0.43701993, 0.07659573)});
	auto trans = RectApprox::Build(left, right, bottom, top);
	add_check(fabs(trans->module()-6.89702)<1e-5, "module");

	Point p(0.519522452, 0.08);
	Point torect = trans->MapToRectangle(p);
	Point topoly = trans->MapToPolygon(torect);
	Point diff = p-topoly;
	add_check(ISZERO(diff.x) && ISZERO(diff.y), "map bnd point: forward, backward");

	p = Point(0.47, 0.02);
	torect = trans->MapToRectangle(p);
	topoly = trans->MapToPolygon(torect);
	diff = p-topoly;
	add_check(ISZERO(diff.x) && ISZERO(diff.y), "map internal point: forward, backward");

}

void test03(){
	using namespace HMMath::Conformal::Impl::DSCPack;
	std::cout<<"DSCPack conformal mapping"<<std::endl; 

	auto bot1 = HMCont2D::Constructor::Circle(8, 2, Point(0,0.1));
	auto top1 = HMCont2D::Constructor::Circle(10, 4, Point(0,0.1));
	auto trans1 = ToAnnulus::Build(top1, bot1);

	auto p10 = vector<Point> {Point(2.5, 1.5)};
	auto p11 = trans1->MapToAnnulus(p10);
	auto p12 = trans1->MapToOriginal(p11);
	add_check(
		ISEQ(trans1->PhiInner(0),trans1->PhiOuter(0)) &&
		ISEQ(trans1->PhiInner(4),trans1->PhiOuter(5)) &&
		p10[0] == p12[0],
		"Uniform points distribution"
	);
}

void test04(){
	std::cout<<"Conformal mapping to rectangle: FEM vs SCPACK"<<std::endl;
	int sz1 = 13;
	auto r1 = HMCont2D::Constructor::Circle(sz1, 10, Point(0,0));
	auto inp1 = HMMath::Conformal::Rect::FactoryInput(r1, {0, sz1/4, sz1/2, 3*sz1/4});
	Point po(7.2, 1.1);

	HMMath::Conformal::Options opt;
	//fem 1
	opt.fem_segment_partition = 12;
	auto trans1 = HMMath::Conformal::Impl::ConfFem::ToRect::Build(
			std::get<0>(inp1), sz1/4, sz1/2, 3*sz1/4, opt);
	Point p1  = trans1->MapToRectangle1(po);
	Point p11 = trans1->MapToPolygon1(p1);

	//fem 2
	opt.fem_segment_partition = 11;
	auto trans2 = HMMath::Conformal::Impl::ConfFem::ToRect::Build(
			std::get<0>(inp1), sz1/4, sz1/2, 3*sz1/4, opt);
	Point p2  = trans2->MapToRectangle1(po);
	Point p22 = trans2->MapToPolygon1(p2);

	//scpack
	auto trans3 = HMMath::Conformal::Impl::SCPack::ToRect::Build(
			std::get<0>(inp1), sz1/4, sz1/2, 3*sz1/4);
	Point p3 =  trans3->MapToRectangle1(po);
	Point p33 = trans3->MapToPolygon1(p3);
			

	add_check(fabs(trans1->module() - 1.05)<0.02, "fem, even partition: module");
	add_check(fabs(trans2->module() - 1.05)<0.02, "fem, odd partition: module");
	add_check(fabs(trans3->module() - 1.05)<0.02, "scpack: module");

	add_check(po == p11, "fem, even: forward->backward");
	add_check(po == p22, "fem, odd: forward->backward");
	add_check(po == p33, "scpack: forward->backward");

	add_check(Point::dist(p1, p3)<0.01, "fem, even vs scpack: forward");
	add_check(Point::dist(p2, p3)<0.01, "fem, odd vs scpack: forward");
}

void test05(){
	std::cout<<"fem in square"<<std::endl;
	Point p1{0.0, 0.0}, p2{2.0, 0.0}, p3{2.0, 1.0}, p4{0.0, 1.0};

	//-> Dfdn(bottom), Dfdx(area), DfDy(area)
	auto dotest = [&](shared_ptr<HMFem::Grid43> g)->std::array<double, 3>{
		auto cont = GGeom::Info::Contour1(*g);
		GridPoint* gp1 = static_cast<GridPoint*>(HMCont2D::ECollection::FindClosestNode(cont, p1));
		GridPoint* gp2 = static_cast<GridPoint*>(HMCont2D::ECollection::FindClosestNode(cont, p2));
		GridPoint* gp3 = static_cast<GridPoint*>(HMCont2D::ECollection::FindClosestNode(cont, p3));
		GridPoint* gp4 = static_cast<GridPoint*>(HMCont2D::ECollection::FindClosestNode(cont, p4));
		auto contbot = HMCont2D::Contour::Assemble(cont, gp1, gp2);
		auto conttop = HMCont2D::Contour::Assemble(cont, gp3, gp4);

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

void test06(){
	std::cout<<"Conformal mapping to annulus: FEM vs DSCPACK"<<std::endl;

	auto top = HMCont2D::Constructor::Circle(12, 1.5, Point(10,0.1));
	auto bot = HMCont2D::Constructor::Circle(5, 0.4, Point(10.5,0.2));
	vector<Point> topp, botp;
	for (auto p: top.ordered_points()) topp.push_back(*p); topp.pop_back();
	for (auto p: bot.ordered_points()) botp.push_back(*p); botp.pop_back();

	Point tp(11.12415, 0.40282);

	//dscpack
	auto trans1 = HMMath::Conformal::Impl::DSCPack::ToAnnulus::Build(topp, botp);
	Point tp1 = trans1->MapToAnnulus1(tp);
	Point tp11 = trans1->MapToOriginal1(tp1);
	add_check(Point::dist(tp, tp11)<1e-2, "dscpack forward->backward");

	//fem
	HMMath::Conformal::Options opt;
	opt.fem_segment_partition = 5;
	auto trans2 = HMMath::Conformal::Impl::ConfFem::ToAnnulus::Build(
			topp, botp, opt);
	Point tp2 = trans2->MapToAnnulus1(tp);
	Point tp22 = trans2->MapToOriginal1(tp2);
	add_check(Point::dist(tp, tp22)<1e-2, "fem forward->backward");

	add_check(fabs(trans1->module() - trans2->module())<1e-2, "modulus");
	add_check([&]()->bool{
		double a1 = ToAngle(trans1->PhiInner(0) - trans1->PhiInner(2));
		double a2 = ToAngle(trans2->PhiInner(0) - trans2->PhiInner(2));
		if (fabs(a1-a2)>1e-1) return false;
		return true;
	}(), "Inner contour angles");
		
	add_check([&]()->bool{
		double a1 = ToAngle(trans1->PhiInner(0) - trans1->PhiOuter(3));
		double a2 = ToAngle(trans2->PhiInner(0) - trans2->PhiOuter(3));
		if (fabs(a1-a2)>1e-1) return false;
		return true;
	}(), "Outer contour angles");
}

int main(){
	test01();
	test02();
	test03();
	test04();
	test05();
	test06();


	if (FAILED_CHECKS ==1){
		std::cout<<FAILED_CHECKS<<" test failed <<<<<<<<<<<<<<<<<<<"<<std::endl;
	} else if (FAILED_CHECKS > 1) {
		std::cout<<FAILED_CHECKS<<" tests failed <<<<<<<<<<<<<<<<<<<"<<std::endl;
	} else {
		std::cout<<"All tests passed"<<std::endl;
	}
	std::cout<<"DONE"<<std::endl;
}
