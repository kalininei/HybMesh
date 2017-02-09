#include "hmblay.hpp"
#include "hmconformal.hpp"
#include "hmtesting.hpp"
#include "constructor.hpp"
#include "cont_assembler.hpp"
#include "infogrid.hpp"
#include "tree.hpp"
#include "export2d_vtk.hpp"
#include "buildgrid.hpp"
#include "modgrid.hpp"
#include "unite_grids.hpp"
#include "finder2d.hpp"

using HMTesting::add_check;

double maxskew(const HM2D::GridData& grid){
	auto s = HM2D::Grid::Skewness(grid);
	return *max_element(s.begin(), s.end());
}

void test01(){
	std::cout<<"01. Boundary layer grids from circle"<<std::endl;
	std::string cn;
	// 1. Build Edges Collection 
	auto col = HM2D::Contour::Constructor::Circle(8, 2.0, Point(0, 0));
	// 2. set inp
	HMBlay::Input inp;
	inp.edges = &col;
	inp.direction = HMBlay::DirectionFromString("OUTER");
	inp.bnd_step_method = HMBlay::MethFromString("NO");
	inp.partition = {0.0, 0.2, 0.4, 1.5};
	inp.bnd_step = 0.1;
	inp.start=inp.end=Point(0,0);
	inp.acute_angle = inp.right_angle = 0.0;
	inp.straight_angle = 300;
	inp.start = inp.end = Point(0,0);
	//3. Calculations
	{
		cn = "Outer full circle";
		HM2D::GridData Ans1 = HMBlay::BuildBLayerGrid({inp});
		add_check(Ans1.vvert.size() == 32 &&
			  Ans1.vcells.size() == 24 &&
			  fabs(HM2D::Grid::Area(Ans1) - 23.35)<0.1, cn);
	}
	{
		cn = "Inner full circle";
		inp.direction = HMBlay::DirectionFromString("INNER");
		HM2D::GridData Ans2 = HMBlay::BuildBLayerGrid({inp});
		add_check(Ans2.vvert.size() == 32 &&
			  Ans2.vcells.size() == 24 &&
			  fabs(HM2D::Grid::Area(Ans2) - 10.903)<0.1, cn);
	}
	{
		cn = "long edges";
		auto col2 = HM2D::Contour::Constructor::Circle(16, 3.0, Point(0, 0));
		HMBlay::Input inp2;
		inp2.edges = &col2;
		inp2.direction = HMBlay::DirectionFromString("OUTER");
		inp2.bnd_step_method = HMBlay::MethFromString("KEEP_SHAPE");
		inp2.partition = {0.0, 0.01, 0.02, 0.03, 0.04};
		inp2.bnd_step = 0.01;
		inp2.start=inp2.end=Point(0,0);
		HM2D::GridData Ans3 = HMBlay::BuildBLayerGrid({inp2});
		add_check( [&]()->bool{
			//points should lie strictly on the source contour
			auto gcont = HM2D::Contour::Tree::GridBoundary(Ans3);
			HM2D::EdgeData& inner = gcont.roots()[0]->children[0].lock()->contour;
			for (auto pt: HM2D::AllVertices(inner)){
				double dist = std::get<4>(HM2D::Contour::CoordAt(col2, *pt));
				if (!ISZERO(dist)) {return false;}
			}
			return true;
		}(), cn);
	}
}

void test02(){
	std::cout<<"02. Smoothing normals"<<std::endl;
	
	//basic inp
	HMBlay::Input inp;
	inp.direction = HMBlay::DirectionFromString("OUTER");
	inp.bnd_step_method = HMBlay::MethFromString("KEEP_SHAPE");
	inp.partition = {0.0, 0.2, 0.4, 1.5};
	inp.bnd_step = 0.05;
	inp.start=inp.end=Point(0,0);
	inp.acute_angle = inp.right_angle = 0.0;
	inp.straight_angle = 300;
	inp.start = inp.end = Point(0,0);
	auto col1 = HM2D::Contour::Constructor::Circle(8, 1.3, Point(2, 2));
	inp.edges = &col1;
	{
		std::string cn("8-side polygon with outer layer");
		HM2D::GridData Ans1 = HMBlay::BuildBLayerGrid({inp});
		vector<double> area1;
		for(int i = 2*Ans1.vcells.size()/3+1; i<Ans1.vcells.size(); ++i)
			area1.push_back(HM2D::Contour::Area(Ans1.vcells[i]->edges));
		//maximum outer element is no more then two times bigger then minimum
		double mina1 = *(std::min_element(area1.begin(), area1.end()));
		double maxa1 = *(std::max_element(area1.begin(), area1.end()));
		add_check(mina1>0 && maxa1/mina1<2.4, cn);
	}
	{
		std::string cn = std::string("8-side polygon with inner layer");
		inp.direction = HMBlay::DirectionFromString("INNER");
		inp.partition = {0.0, 0.02, 0.08, 0.2};
		HM2D::GridData Ans2 = HMBlay::BuildBLayerGrid({inp});
		vector<double> area2;
		for(int i = 0; i<Ans2.vcells.size()/3; ++i)
			area2.push_back(HM2D::Contour::Area(Ans2.vcells[i]->edges));
		//maximum outer element is no more than two times larger than minimum
		double mina2 = *(std::min_element(area2.begin(), area2.end()));
		double maxa2 = *(std::max_element(area2.begin(), area2.end()));
		add_check(mina2>0 && maxa2/mina2<2.4, cn);
	}
};


void test03(){
	std::cout<<"03. Open contour boundary layer"<<std::endl;
	
	//basic inp
	HMBlay::Input inp;
	inp.bnd_step_method = HMBlay::MethFromString("KEEP_ALL");
	inp.partition = {0.0, 0.2, 0.4, 1.5};
	inp.bnd_step = 0.3;
	inp.start=inp.end=Point(0,0);
	inp.acute_angle = inp.right_angle = 0.0;
	inp.straight_angle = 300;

	std::string cn("5-points polyline. Left layer");
	auto col1 = HM2D::Contour::Constructor::FromPoints(
			{0,0, 1,1, 2,1.9, 4,3.5, 7,6.1});
	inp.direction = HMBlay::DirectionFromString("LEFT");
	inp.start = Point(0,0);
	inp.end = Point(7,6.1);
	inp.edges = &col1;
	HM2D::GridData Ans1 = HMBlay::BuildBLayerGrid({inp});
	add_check([&](){
		for (int i=0; i<Ans1.vcells.size(); ++i) if (HM2D::Contour::Area(Ans1.vcells[i]->edges)<0) return false;
		for (int i=0; i<Ans1.vvert.size(); ++i) if (Ans1.vvert[i]->y<0.0) return false;
		return true;
	}(), cn);

	cn = std::string("5-points polyline. Right layer");
	inp.direction = HMBlay::DirectionFromString("RIGHT");
	HM2D::GridData Ans2 = HMBlay::BuildBLayerGrid({inp});
	add_check([&](){
		for (int i=0; i<Ans2.vcells.size(); ++i) if (HM2D::Contour::Area(Ans2.vcells[i]->edges)<0) return false;
		for (int i=0; i<Ans2.vvert.size(); ++i) if (Ans2.vvert[i]->y>6.1) return false;
		return true;
	}(), cn);

	cn = std::string("5-points polyline. Right layer. Opposite direction");
	std::swap(inp.start, inp.end);
	HM2D::GridData Ans3 = HMBlay::BuildBLayerGrid({inp});
	HM2D::Export::GridVTK(Ans3, "g2.vtk");
	add_check([&](){
		for (int i=0; i<Ans3.vcells.size(); ++i) if (HM2D::Contour::Area(Ans3.vcells[i]->edges)<0) return false;
		for (int i=0; i<Ans3.vvert.size(); ++i) if (Ans3.vvert[i]->y<0.0) return false;
		return true;
	}(), cn);
}

void test04(){
	std::cout<<"04. Different partitions"<<std::endl;

	auto col1 = HM2D::Contour::Constructor::Circle(8, 10, Point(0, 0));
	//basic inp1
	HMBlay::Input inp1;
	inp1.bnd_step_method = HMBlay::MethFromString("KEEP_ALL");
	inp1.acute_angle = inp1.right_angle = 0.0;
	inp1.straight_angle = 300;
	inp1.direction = HMBlay::DirectionFromString("OUTER");
	inp1.edges = &col1;

	std::string cn = std::string("mesh half of closed contour");
	inp1.partition = {0, 0.1, 0.2, 0.3, 0.4};
	inp1.start = Point(10, 0);
	inp1.end = Point(-10, 0);
	inp1.bnd_step = 0.3;
	HM2D::GridData Ans1 = HMBlay::BuildBLayerGrid({inp1});
	add_check([&](){
		for (int i=0; i<Ans1.vvert.size(); ++i) if (Ans1.vvert[i]->y<-1e-3) return false;
		bool hasp1=false, hasp2=false;
		for (int i=0; i<Ans1.vvert.size(); ++i){
			if (Point::dist(*Ans1.vvert[i], Point(10.4, 0))<1e-2) hasp1 = true;
			if (Point::dist(*Ans1.vvert[i],Point(-10.4, 0))<1e-2) hasp2 = true;
		}
		return true;
	}(), cn);

	cn = std::string("Two layers with same partition depth count: outer");
	HMBlay::Input inp2(inp1);
	inp2.partition = {0, 0.1, 0.2, 0.3, 0.8};
	std::swap(inp2.start, inp2.end);
	inp2.bnd_step = 0.1;
	HM2D::GridData Ans2 = HMBlay::BuildBLayerGrid({inp1, inp2});
	add_check([&](){
		for (int i=0; i<Ans2.vcells.size(); ++i) if (HM2D::Contour::Area(Ans2.vcells[i]->edges)<0) return false;
		if (Ans2.vcells.size()/4 != Ans2.vvert.size()/5) return false;
		for (int i=0; i<Ans2.vvert.size(); ++i) if (Ans2.vvert[i]->y<-10-0.8-1e-12) return false;
		for (int i=0; i<Ans2.vvert.size(); ++i) if (Ans2.vvert[i]->y>10+0.4+1e-12) return false;
		bool hasp1=false, hasp2=false;
		for (int i=0; i<Ans2.vvert.size(); ++i){
			if (Point::dist(*Ans2.vvert[i], Point(10.8, 0))<5e-2) hasp1 = true;
			if (Point::dist(*Ans2.vvert[i],Point(-10.8, 0))<5e-2) hasp2 = true;
		}
		if (!hasp1 || !hasp2) return false;
		return true;
	}(), cn);
	HM2D::Export::GridVTK(Ans2, "g1.vtk");

	cn = std::string("Two layers with same partition depth count: inner");
	inp1.direction = inp2.direction = HMBlay::DirectionFromString("INNER");
	HM2D::GridData Ans3 = HMBlay::BuildBLayerGrid({inp1, inp2});
	add_check([&](){
		for (int i=0; i<Ans3.vcells.size(); ++i) if (HM2D::Contour::Area(Ans3.vcells[i]->edges)<0) return false;
		if (Ans3.vcells.size()/4 != Ans3.vvert.size()/5) return false;
		for (int i=0; i<Ans3.vvert.size(); ++i){
			const HM2D::Vertex* p=Ans3.vvert[i].get();
			if (p->y>10+1e-3 || p->y<-10-1e-3) return false;
			if (ISEQ(p->x, 0) && (p->y>-9.13 && p->y<9.45)) return false;
		}
		return true;
	}(), cn);

	cn = std::string("Two layers with different partition depth count");
	inp1.direction = inp2.direction = HMBlay::DirectionFromString("OUTER");
	inp2.partition.push_back(1.0); inp2.partition.push_back(1.2);
	HM2D::GridData Ans4 = HMBlay::BuildBLayerGrid({inp1, inp2});
	add_check([&](){
		for (int i=0; i<Ans4.vcells.size(); ++i) if (HM2D::Contour::Area(Ans4.vcells[i]->edges)<0) return false;
		int cp = 0;
		for (int i=0; i<Ans4.vvert.size(); ++i){
			if (fabs(Ans4.vvert[i]->y)<0.05){
				if (Ans4.vvert[i]->x >  10.81) ++cp;
				if (Ans4.vvert[i]->x < -10.81) ++cp;
			}
		}
		if (cp!=4) return false;
		return true;
	}(), cn);
}


void test05(){
	std::cout<<"05. Non-right boundary angle at the ends"<<std::endl;
	auto col1 = HM2D::Contour::Constructor::FromPoints(
		{0,0, 5,0, 6,4, 2,4}, true);

	HMBlay::Input inp1;
	inp1.bnd_step_method = HMBlay::MethFromString("KEEP_ALL");
	inp1.direction = HMBlay::DirectionFromString("INNER");
	inp1.edges = &col1;
	inp1.partition = {0, 0.1, 0.2, 0.3, 0.4};
	inp1.start = Point(0, 0);
	inp1.end = Point(5, 0);
	inp1.bnd_step = 0.1;

	std::string cn = std::string("mesh half of closed contour");
	HM2D::GridData Ans1 = HMBlay::BuildBLayerGrid({inp1});
	add_check([&](){
		if (fabs(HM2D::Grid::Area(Ans1)-1.93)>0.01) return false;
		for (int i=0; i<Ans1.vcells.size(); ++i) if (HM2D::Contour::Area(Ans1.vcells[i]->edges)<0) return false;
		return true;
	}(), cn);
}


void test06(){
	std::cout<<"06. 1 segment sources"<<std::endl;
	auto col1 = HM2D::Contour::Constructor::FromPoints(
			{0,0, 2,0});

	HMBlay::Input inp1;
	inp1.bnd_step_method = HMBlay::MethFromString("NO");
	inp1.acute_angle = inp1.right_angle = 0.0;
	inp1.straight_angle = 300;
	inp1.partition = {0, 0.1};
	inp1.bnd_step = 0.1;

	std::string cn = std::string("left from positive");
	inp1.edges = &col1;
	inp1.direction = HMBlay::DirectionFromString("LEFT");
	inp1.start = Point(0, 0);
	inp1.end = Point(2, 0);
	HM2D::GridData Ans1 = HMBlay::BuildBLayerGrid({inp1});
	add_check([&](){
		if (Ans1.vcells.size() != 1) return false;
		for (int i=0; i<Ans1.vvert.size(); ++i) if (Ans1.vvert[i]->y<-1e-12) return false;
		return true;
	}(), cn);

	cn = std::string("right from positive");
	inp1.direction = HMBlay::DirectionFromString("RIGHT");
	HM2D::GridData Ans2 = HMBlay::BuildBLayerGrid({inp1});
	add_check([&](){
		if (Ans2.vcells.size() != 1) return false;
		for (int i=0; i<Ans2.vvert.size(); ++i) if (Ans2.vvert[i]->y>1e-12) return false;
		return true;
	}(), cn);

	cn = std::string("left from negative");
	std::swap(inp1.start, inp1.end);
	inp1.direction = HMBlay::DirectionFromString("LEFT");
	HM2D::GridData Ans3 = HMBlay::BuildBLayerGrid({inp1});
	add_check([&](){
		if (Ans3.vcells.size() != 1) return false;
		for (int i=0; i<Ans3.vvert.size(); ++i) if (Ans3.vvert[i]->y>1e-12) return false;
		return true;
	}(), cn);

	cn = std::string("right from negative");
	inp1.direction = HMBlay::DirectionFromString("RIGHT");
	HM2D::GridData Ans4 = HMBlay::BuildBLayerGrid({inp1});
	add_check([&](){
		if (Ans4.vcells.size() != 1) return false;
		for (int i=0; i<Ans4.vvert.size(); ++i) if (Ans4.vvert[i]->y<-1e-12) return false;
		return true;
	}(), cn);

	cn = std::string("right from negative. fine partition");
	inp1.bnd_step_method = HMBlay::MethFromString("IGNORE_ALL");
	HM2D::GridData Ans5 = HMBlay::BuildBLayerGrid({inp1});
	add_check([&](){
		if (Ans5.vcells.size() != 20) return false;
		for (int i=0; i<Ans5.vvert.size(); ++i) if (Ans5.vvert[i]->y<-1e-12) return false;
		return true;
	}(), cn);
}


void test07(){
	std::cout<<"07. Corner angle algorithm basic test"<<std::endl;
	auto c1 = HM2D::Contour::Constructor::Circle(36, 3, Point(0, 0));
	auto vc1 = HM2D::AllVertices(c1);
	HM2D::VertexData apoints;
	apoints.emplace_back(new HM2D::Vertex(-3, -7));
	apoints.emplace_back(new HM2D::Vertex(4, -7));
	Point *p0, *p1;
	p0 = vc1[std::get<0>(HM2D::Finder::ClosestPoint(vc1, Point(3,0)))].get();
	p1 = vc1[std::get<0>(HM2D::Finder::ClosestPoint(vc1, Point(-3,0)))].get();
	HM2D::EdgeData con1 = HM2D::Contour::Assembler::ShrinkContour(c1, p0, p1);
	con1.emplace_back(new HM2D::Edge(HM2D::Contour::Last(con1), apoints[0]));
	con1.emplace_back(new HM2D::Edge(apoints[0], apoints[1]));
	con1.emplace_back(new HM2D::Edge(apoints[1], HM2D::Contour::First(con1)));
	
	HMBlay::Input inp1;
	inp1.bnd_step_method = HMBlay::MethFromString("KEEP_SHAPE");
	inp1.acute_angle = 0;
	inp1.right_angle = 120;
	inp1.straight_angle = 235;
	inp1.direction = HMBlay::DirectionFromString("INNER");
	inp1.partition = {0, 0.1, 0.2, 0.3, 0.4};
	inp1.start = inp1.end = Point(0, 0);
	inp1.bnd_step = 0.3;
	inp1.edges = &con1;
	HM2D::GridData Ans1 = HMBlay::BuildBLayerGrid({inp1});
	add_check(Ans1.vvert.size() == 545 && Ans1.vcells.size() == 436,
			"domain with a half-cirlce");

}

void test08(){
	std::cout<<"08 Single angle connection"<<std::endl;
	auto c1 = HM2D::Contour::Constructor::Circle(36, 3, Point(0, 0));
	auto vc1 = HM2D::AllVertices(c1);
	HM2D::VertexData apoints;
	apoints.emplace_back(new HM2D::Vertex(3, -3));
	Point *p0, *p1;
	p0 = vc1[std::get<0>(HM2D::Finder::ClosestPoint(vc1, Point(3,0)))].get();
	p1 = vc1[std::get<0>(HM2D::Finder::ClosestPoint(vc1, Point(0,-3)))].get();
	HM2D::EdgeData con1 = HM2D::Contour::Assembler::ShrinkContour(c1, p0, p1);
	con1.emplace_back(new HM2D::Edge(HM2D::Contour::Last(con1), apoints[0]));
	con1.emplace_back(new HM2D::Edge(apoints[0], HM2D::Contour::First(con1)));
	HM2D::Export::ContourVTK(con1, "c1.vtk");

	HMBlay::Input inp1;
	inp1.bnd_step_method = HMBlay::MethFromString("KEEP_SHAPE");
	inp1.direction = HMBlay::DirectionFromString("INNER");
	inp1.partition = {0, 0.02, 0.05, 0.1, 0.2};
	inp1.start = inp1.end = Point(0, 0);
	inp1.bnd_step = 0.2;
	inp1.edges = &con1;

	HM2D::GridData Ans1 = HMBlay::BuildBLayerGrid({inp1});
	auto sz1 = HM2D::Grid::CellAreas(Ans1);
	double minsz1 = *std::min_element(sz1.begin(), sz1.end());
	double maxsz1 = *std::max_element(sz1.begin(), sz1.end());
	add_check(Ans1.vvert.size() == 565 && Ans1.vcells.size() == 452 &&
	          fabs(maxsz1 - 0.02)<1e-6 && fabs(minsz1 - 0.0004)<1e-6, "Mesh inner");

	inp1.direction = HMBlay::DirectionFromString("OUTER");
	HM2D::GridData Ans2 = HMBlay::BuildBLayerGrid({inp1});
	auto sz2 = HM2D::Grid::CellAreas(Ans2);
	double minsz2 = *std::min_element(sz2.begin(), sz2.end());
	double maxsz2 = *std::max_element(sz2.begin(), sz2.end());
	add_check(Ans2.vvert.size() == 575 && Ans2.vcells.size() == 460 &&
	          fabs(maxsz2 - 0.0209396)<1e-4 && fabs(minsz2 - 0.0004)<1e-6, "Mesh outer");
	HM2D::Export::GridVTK(Ans1, "g1.vtk");
	HM2D::Export::GridVTK(Ans2, "g2.vtk");
}

void test09(){
	std::cout<<"09 Throw on crossing intervals"<<std::endl;
	//TODO
	//HMFem::Mat m;
	//for (int i=0; i<100; ++i) m.set(i,i,1);
	//for (int i=0; i<99; ++i) m.set(i, i+1, 0.3);
	//auto slv = HMFem::MatSolve::Factory(m);
}

void test10(){
	std::cout<<"10. right + re-entrant angle"<<std::endl;
	auto c = HM2D::Contour::Constructor::FromPoints(
		{0,0, 2,0, 2.3,-1,  5,-1});
	HMBlay::Input inp1;
	inp1.bnd_step_method = HMBlay::MethFromString("KEEP_SHAPE");
	inp1.acute_angle = 0;
	inp1.right_angle = 150;
	inp1.straight_angle = 200;
	inp1.direction = HMBlay::DirectionFromString("INNER");
	inp1.partition = {0, 0.05, 0.2, 0.3, 0.4};
	inp1.start = Point(0, 0);
	inp1.end = Point(5, -1);
	inp1.bnd_step = 0.1;
	inp1.edges = &c;
	inp1.force_conformal = false;
	HM2D::GridData Ans1 = HMBlay::BuildBLayerGrid({inp1});
	add_check(Ans1.vcells.size() == 232 && Ans1.vvert.size() == 295, "Mesh");
};

void test11(){
	std::cout<<"11. Acute angle"<<std::endl;
	auto c = HM2D::Contour::Constructor::FromPoints(
			{-10,0, 5,0, 0,1.5});
	HMBlay::Input inp1;
	inp1.bnd_step_method = HMBlay::MethFromString("KEEP_SHAPE");
	inp1.direction = HMBlay::DirectionFromString("INNER");
	inp1.edges = &c;
	inp1.bnd_step = 0.1;
	inp1.start = Point(0,0);
	inp1.end = Point(0,1.5);
	inp1.partition = {0, 0.01, 0.05, 0.1, 0.15, 0.2};

	{
		HM2D::GridData Ans1 = HMBlay::BuildBLayerGrid({inp1});
		int badskew1=0;
		for (double s: HM2D::Grid::Skewness(Ans1)) {if (s>0.7) ++badskew1;}
		add_check(Ans1.vcells.size() > 550 && Ans1.vvert.size() > 590 && badskew1==1, "Mesh1");
		HM2D::Export::GridVTK(Ans1, "g1.vtk");
	}
	{
		inp1.partition = {0, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2};
		HM2D::GridData Ans2 = HMBlay::BuildBLayerGrid({inp1});
		int badskew2=0;
		for (double s: HM2D::Grid::Skewness(Ans2)) {if (s>0.7) ++badskew2;}
		add_check(Ans2.vcells.size() > 700 && Ans2.vvert.size() > 720 && badskew2==1, "Mesh2");
	}

};

void test12(){
	std::cout<<"12. Different kinds of reflex angles"<<std::endl;
	auto c = HM2D::Contour::Constructor::FromPoints(
			{0,0, 1,0.3, 2,0, 1.8,-3, 1.5,-1});

	HMBlay::Input inp1;
	inp1.bnd_step_method = HMBlay::MethFromString("KEEP_SHAPE");
	inp1.direction = HMBlay::DirectionFromString("LEFT");
	inp1.edges = &c;
	inp1.bnd_step = 0.05;
	inp1.start = Point(0,0);
	inp1.end = Point(1.5,-1);
	inp1.partition = {0, 0.01, 0.02, 0.04, 0.07, 0.1};
	inp1.acute_angle = 30;
	inp1.right_angle = 100;
	inp1.straight_angle = 240;
	inp1.reentrant_angle = 300;

	HM2D::GridData Ans1 = HMBlay::BuildBLayerGrid({inp1});
	add_check(Ans1.vvert.size() == 938 && Ans1.vcells.size() == 785, "Mesh");
}

void test13(){
	std::cout<<"13. Cylinder in channel meshing"<<std::endl;
	//g1
	auto g1 = HM2D::Grid::Constructor::RectGrid(Point(0,0), Point(30, 10), 90, 30);
	//g2
	HMBlay::Input inp1, inp2;
	auto c = HM2D::Contour::Constructor::FromPoints({0, 0, 30, 0});
	inp1.bnd_step_method = HMBlay::MethFromString("KEEP_SHAPE");
	inp1.direction = HMBlay::DirectionFromString("LEFT");
	inp1.edges = &c;
	inp1.bnd_step = 0.2;
	inp1.start = Point(0,0);
	inp1.end = Point(30, 0);
	inp1.partition = {0, 0.03, 0.07, 0.12, 0.2, 0.3, 0.45, 0.7};
	HM2D::GridData g2_1 = HMBlay::BuildBLayerGrid({inp1});
	c = HM2D::Contour::Constructor::FromPoints({30, 10, 0, 10});
	inp1.edges = &c;
	inp1.start = Point(30,10); inp1.end = Point(0, 10);
	HM2D::GridData g2_2 = HMBlay::BuildBLayerGrid({inp1});
	HM2D::GridData g2 = g2_1;
	HM2D::Grid::Algos::ShallowAdd(g2_2, g2);
	//g3
	c = HM2D::Contour::Constructor::Circle(64, 0.5, Point(10, 5));
	inp1.direction = HMBlay::DirectionFromString("RIGHT");
	inp1.edges = &c;
	inp1.bnd_step = 0.05;
	inp1.partition = {0, 0.03, 0.07, 0.12, 0.18, 0.24};
	inp1.start = Point(10.6, 7);
	inp1.end = Point(10.6, 3);
	inp2 = inp1;
	inp2.bnd_step = 0.02;
	inp2.start = Point(10.6, 3);
	inp2.end = Point(10.6, 7);
	inp2.partition = {0, 0.03, 0.07, 0.12, 0.18, 0.24, 0.3, 0.36, 0.42, 0.48, 0.54, 0.6};
	HM2D::GridData g3 = HMBlay::BuildBLayerGrid({inp1, inp2});
	//g4
	auto conf = HMMap::Conformal::Rect::Factory(
		{Point(10.38, 4.07), Point(15.33, 3), Point(15.33, 7), Point(10.38, 5.93)},
		{0,1,2,3});
	auto g4 = HM2D::Grid::Constructor::RectGrid(Point(0, 0), Point(conf->module(), 1), 30, 50);
	for (auto& v: g4.vvert){
		Point a = conf->MapToPolygon1(*v);
		v->set(a);
	}
	//g5
	auto g5 = HM2D::Grid::Constructor::RectGrid(Point(14.5, 2.5), Point(30, 7.5), 92, 30);

	//gcommon
	HM2D::Grid::Algos::OptUnite opt;
	opt.buffer_size = 0.3;
	auto gr1 = HM2D::Grid::Algos::UniteGrids.ToCout(g1, g5, opt);
	opt.buffer_size = 0.1;
	auto gr2 = HM2D::Grid::Algos::UniteGrids.ToCout(gr1, g4, opt);
	auto gr3 = HM2D::Grid::Algos::UniteGrids.ToCout(gr2, g2, opt);
	opt.buffer_size = 0.2;
	opt.empty_holes = true;
	auto gr4 = HM2D::Grid::Algos::UniteGrids.ToCout(gr3, g3, opt);
	add_check(maxskew(gr4)<0.8, "skewness check");
}

HM2D::GridData grid_minus_cont(HM2D::GridData& g, const HM2D::EdgeData& cont){
	vector<Point> pp;
	vector<int> eds;
	for (auto p: HM2D::Contour::OrderedPoints(cont)) pp.push_back(*p);
	pp.resize(pp.size()-1);
	for (int i=0; i<pp.size(); ++i){
		eds.push_back(i);
		eds.push_back(i+1);
	}
	eds.back() = 0;
	auto tree = HM2D::Contour::Tree::Assemble(cont);
	return HM2D::Grid::Algos::SubstractArea.ToCout(g, tree, false);
}

void test14(){
	std::cout<<"14. Four angle types"<<std::endl;

	//right angle
	auto fx = [](double t){ return 0.01*sin(2.0*M_PI*3.0*t); };
	auto fy = [](double t){ return 0.013*(-sin(2.0*M_PI*5.0*t)); };
	double h=0.01;
	double t=1-h;
	HM2D::VertexData pc1;
	while (t>h/2.0){
		pc1.emplace_back(new HM2D::Vertex(fx(t)+0.08*t, t));
		t-=h;
	}
	while (t<1+h/2.0){
		pc1.emplace_back(new HM2D::Vertex(t, fy(t)));
		t+=h;
	}
	pc1.emplace_back(new HM2D::Vertex(1,1));
	auto c1 = HM2D::Contour::Assembler::Contour1(pc1, true); 
	HMBlay::Input inp1;
	inp1.bnd_step_method = HMBlay::MethFromString("IGNORE_ALL");
	inp1.direction = HMBlay::DirectionFromString("LEFT");
	inp1.edges = &c1;
	inp1.bnd_step = 0.03;
	inp1.force_conformal = true;
	inp1.partition = {0, 0.01, 0.02, 0.03, 0.04, 0.053, 0.07, 0.095, 0.12};

	HM2D::Export::ContourVTK(c1, "ecol.vtk");

	HM2D::GridData g1 = HMBlay::BuildBLayerGrid({inp1});
	HM2D::GridData basgrid1 = HM2D::Grid::Constructor::RectGrid(Point(-0.1, -0.1), Point(1.1, 1.1), 30, 30);
	HM2D::GridData b1 = grid_minus_cont(basgrid1, HM2D::Contour::Assembler::GridBoundary1(g1));
	HM2D::Grid::Algos::OptUnite opt(0.03);
	HM2D::GridData gr1 = HM2D::Grid::Algos::UniteGrids.ToCout(b1, g1, opt); 
	HM2D::Grid::Algos::SimplifyBoundary(gr1, 45);
	add_check(fabs(HM2D::Grid::Area(gr1) - 0.955398)<1e-5 &&
		  maxskew(gr1)<0.7, "Square with curved boundary");

	//reentrant
	vector<Point> pc2;
	for (auto an = -M_PI/2.0; an<=0.0; an+=M_PI/2.0/20){
		pc2.push_back(Point(cos(an)-1, sin(an)));
	}
	for (auto an = -M_PI/2.0; an<=0.0; an+=M_PI/2.0/20){
		pc2.push_back(Point(cos(an), sin(an)));
	}
	pc2.push_back(Point(4, 0.4));
	pc2.push_back(Point(4, 1.5));
	pc2.push_back(Point(-1, 1.0));
	auto c2 = HM2D::Contour::Constructor::FromPoints(pc2, true);
	inp1.edges=&c2;
	inp1.bnd_step = 0.02;
	inp1.start =inp1.end = Point(0, 0);
	HM2D::GridData g2 = HMBlay::BuildBLayerGrid({inp1});
	HM2D::GridData basgrid2 = HM2D::Grid::Constructor::RectGrid(Point(-2.1, -1.1), Point(2.8, 2.1), 150, 100);
	HM2D::GridData b2 = grid_minus_cont(basgrid2, HM2D::Contour::Assembler::GridBoundary1(g2));
	HM2D::Grid::Algos::OptUnite opt2(0.03);
	HM2D::GridData gr2 = HM2D::Grid::Algos::UniteGrids.ToCout(b2, g2, opt2);
	HM2D::Grid::Algos::SimplifyBoundary(gr2, 45);
	add_check(fabs(HM2D::Grid::Area(gr2) - 6.26758)<1e-5, "Reentrant area with a hole");

	//acute angle
	{
		auto c3 = HM2D::Contour::Constructor::FromPoints({0,0, 5,0, 0,1.5}, true);
		HMBlay::Input inp2;
		inp2.bnd_step_method = HMBlay::MethFromString("KEEP_SHAPE");
		inp2.direction = HMBlay::DirectionFromString("INNER");
		inp2.edges = &c3;
		inp2.bnd_step = 0.07;
		inp2.force_conformal = true;
		inp2.start = inp2.end = Point(0,0);
		inp2.partition = {0, 0.02, 0.05, 0.1, 0.15, 0.2};
		HM2D::GridData g3 = HMBlay::BuildBLayerGrid({inp2});
		HM2D::GridData basgrid3 = HM2D::Grid::Constructor::RectGrid(Point(-0.1, -0.1), Point(4.5, 2), 80 , 40);
		HM2D::Grid::Algos::OptUnite opt3(0.03, true);
		HM2D::GridData r3 = HM2D::Grid::Algos::UniteGrids.ToCout(basgrid3, g3, opt3);
		HM2D::GridData gr3 = grid_minus_cont(r3, c3);
		add_check(fabs(HM2D::Grid::Area(gr3) - 3.75)<1e-5, "Closed tringle with acute angle");
	}
};

void test15(){
	std::cout<<"15. Geometry noise"<<std::endl;

	auto c1 = HM2D::Contour::Constructor::FromPoints({0,0, 0.5, -0.03, 0.55, 0, 0.68, 0, 0.7, 0.02, 0.75, 0, 1, 0, 1.01, -0.02, 1.03, 0, 1.1, 0}); 
	HM2D::Export::ContourVTK(c1, "c15.vtk");
	HMBlay::Input inp1;
	inp1.bnd_step_method = HMBlay::MethFromString("KEEP_SHAPE");
	inp1.direction = HMBlay::DirectionFromString("LEFT");
	inp1.edges = &c1;
	inp1.bnd_step = 0.05;
	inp1.start = Point(0, 0);
	inp1.end = Point(10, 0);
	inp1.partition = {0, 0.01, 0.02, 0.04, 0.08, 0.12, 0.16, 0.2};
	HM2D::GridData g1 = HMBlay::BuildBLayerGrid({inp1});

}

void test16(){
	std::cout<<"16. Impostition of a boundary grid with contour conflicts"<<std::endl;
	//1) cicle grid
	HM2D::GridData gcirc = HM2D::Grid::Constructor::Ring(Point(0, 0), 1, 0.6, 34, 5);

	//2) get boundary from it
	HM2D::Contour::Tree gcont = HM2D::Contour::Tree::GridBoundary(gcirc);
	auto _ae = gcont.alledges();
	//3) build a boundary layer around it
	HMBlay::Input inp1;
	inp1.bnd_step_method = HMBlay::MethFromString("IGNORE_ALL");
	inp1.direction = HMBlay::DirectionFromString("LEFT");
	inp1.edges = &_ae;
	inp1.bnd_step = 0.06;
	inp1.start = Point(1, 0);
	inp1.end = Point(1, 1);
	inp1.partition = {0, 0.01, 0.02, 0.04};

	HM2D::GridData bgrid = HMBlay::BuildBLayerGrid({inp1}); 
	add_check(bgrid.vvert.size() == 60, "IGNORE_ALL boundary partition");

	//simplify boundary tests
	HM2D::GridData bgrid2;
	HM2D::DeepCopy(bgrid, bgrid2);
	HM2D::Grid::Algos::SimplifyBoundary(bgrid2, 180);
	add_check(bgrid2.vvert.size() == 52, "full simplification of bgrid");

	HM2D::GridData bgrid4;
	HM2D::DeepCopy(bgrid, bgrid4);
	HM2D::Grid::Algos::SimplifyBoundary(bgrid4, 45);
	add_check(bgrid4.vvert.size()==56, "simplification with angle=pi/4");

	HM2D::GridData bgrid5;
	HM2D::DeepCopy(bgrid, bgrid5);
	HM2D::Grid::Algos::SimplifyBoundary(bgrid5, 0.0);
	add_check(fabs(HM2D::Grid::Area(bgrid5)-HM2D::Grid::Area(bgrid))<1e-6, "simplification with angle=0");

	//impose
	HM2D::Grid::Algos::OptUnite opt(0.3);
	HM2D::GridData impgrid = HM2D::Grid::Algos::UniteGrids.ToCout(gcirc, bgrid, opt);
	HM2D::GridData impgrid1 = HM2D::Grid::Algos::UniteGrids.ToCout(gcirc, bgrid4, opt);

	//see the result
	double arinit = HM2D::Grid::Area(gcirc);
	add_check(ISEQ(arinit, HM2D::Grid::Area(impgrid)), "imposition with non-simplified bgrid: area");
	add_check(ISEQ(arinit, HM2D::Grid::Area(impgrid1)), "imposition with simplified bgrid: area");
	HM2D::Grid::Algos::SimplifyBoundary(impgrid, 45);
	HM2D::Grid::Algos::SimplifyBoundary(impgrid1, 45);
	add_check((fabs(1.99894 - HM2D::Grid::Area(impgrid))<1e-4), "area1 after simplification");
	add_check(ISEQ(arinit, HM2D::Grid::Area(impgrid1)), "area2 after simplification. Keep all cell non-degenerate");
	
	auto skew1 = HM2D::Grid::Skewness(impgrid);
	auto skew2 = HM2D::Grid::Skewness(impgrid1);
	double maxskew1 = *std::max_element(skew1.begin(), skew1.end());
	double maxskew2 = *std::max_element(skew2.begin(), skew2.end());
	vector<double> badskew1, badskew2;
	std::copy_if(skew1.begin(), skew1.end(), std::back_inserter(badskew1), [](double a){return a>0.8;});
	std::copy_if(skew2.begin(), skew2.end(), std::back_inserter(badskew2), [](double a){return a>0.8;});
	add_check(badskew1.size() == 0 && badskew2.size() < 6, "skewness check");
	HM2D::Export::GridVTK(impgrid, "g1.vtk");
	HM2D::Export::GridVTK(impgrid1, "g2.vtk");
}

void test17(){
	std::cout<<"17. Snapping to source contour"<<std::endl;
	std::vector<Point> inpx, inpy;
	for (int i=1; i<20; ++i) inpx.push_back(Point(0.1*i, 0));
	for (int i=1; i<20; ++i) inpy.push_back(Point(0, 0.1*i));
	double sx = 0.02, sy = 0.02;
	int coef=1;
	for (auto& p: inpx) { p.y += coef * sx; coef *= -1; }
	for (auto& p: inpy) { p.x += coef * sy; coef *= -1; }

	std::vector<Point> inppts;
	inppts.insert(inppts.end(), inpy.rbegin(), inpy.rend());
	inppts.push_back(Point(0, 0));
	inppts.insert(inppts.end(), inpx.begin(), inpx.end());

	HM2D::EdgeData
	cont = HM2D::Contour::Constructor::FromPoints(inppts, false);

	HMBlay::Input inp1;
	inp1.bnd_step_method = HMBlay::MethFromString("NO");
	inp1.direction = HMBlay::DirectionFromString("LEFT");
	inp1.edges = &cont;
	inp1.start = Point(0, 0);
	inp1.end = Point(2, 0);
	inp1.partition = {0, 0.01, 0.02, 0.04, 0.08, 0.14, 0.2};

	HM2D::GridData bgrid1 = HMBlay::BuildBLayerGrid({inp1}); 
	auto bcont1 = HM2D::Contour::Tree::GridBoundary(bgrid1);
	add_check(bgrid1.vvert.size() == 141 &&
		bcont1.whereis(Point(-0.02+1e-3, 0.1)) == INSIDE &&
		bcont1.whereis(Point(1.88046632081399, 0.209839818196609)) == INSIDE,
		"one section with zig-zag bottom and left");
	HM2D::Export::ContourVTK(cont, "t17_1.vtk");

	inp1.start = Point(-0.02, 1.9);
	HM2D::GridData bgrid2 = HMBlay::BuildBLayerGrid({inp1});
	auto bcont2 = HM2D::Contour::Tree::GridBoundary(bgrid2);
	add_check(bgrid2.vvert.size() == 297 &&
		bcont2.whereis(Point(-0.02+1e-3, 0.1)) == INSIDE &&
		bcont2.whereis(Point(0.1, 0.0168026937561729)) == OUTSIDE,
		"zig-zag bottom and left sections with right angle");

	inp1.bnd_step_method = HMBlay::MethFromString("IGNORE_ALL");
	inp1.bnd_step = 0.04;
	HM2D::GridData bgrid3 = HMBlay::BuildBLayerGrid({inp1});
	auto bcont3 = HM2D::Contour::Tree::GridBoundary(bgrid2);
	add_check(bgrid3.vvert.size() == 722 &&
		bcont3.whereis(Point(-0.02+1e-3, 0.1)) == INSIDE &&
		bcont3.whereis(Point(0.1, 0.0168026937561729)) == OUTSIDE && 
		bcont3.whereis(Point(0.0161764309193625, 1.80159049128253)) == OUTSIDE,
		"same with ignore_all option");
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
	test11();
	test12();
	test13();
	test14();
	test16();
	test17();
	
	//UNDONE:
	//test15();
	
	HMTesting::check_final_report();
	std::cout<<"DONE"<<std::endl;
}
