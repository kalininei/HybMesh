#include "hmblay.hpp"
#include "fileproc.h"
#include "hmconformal.hpp"
#include "procgrid.h"

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
	std::cout<<"01. Boundary Layers grid from circle"<<std::endl;
	std::string cn;
	// 1. Build Edges Collection 
	auto col = HMCont2D::Constructor::Circle(8, 2.0, Point(0, 0));
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
	cn = "Outer full circle";
	GridGeom Ans1 = HMBlay::BuildBLayerGrid({inp});
	add_check(Ans1.n_points() == 32 &&
		  Ans1.n_cells() == 24 &&
		  fabs(Ans1.area() - 23.35)<0.1, cn);

	cn = "Inner full circle";
	inp.direction = HMBlay::DirectionFromString("INNER");
	GridGeom Ans2 = HMBlay::BuildBLayerGrid({inp});
	add_check(Ans2.n_points() == 32 &&
		  Ans2.n_cells() == 24 &&
		  fabs(Ans2.area() - 10.903)<0.1, cn);
	save_vtk(Ans1, "t01_1.vtk");
	save_vtk(Ans2, "t01_2.vtk");

	cn = "long edges";
	auto col2 = HMCont2D::Constructor::Circle(16, 3.0, Point(0, 0));
	HMBlay::Input inp2;
	inp2.edges = &col2;
	inp2.direction = HMBlay::DirectionFromString("OUTER");
	inp2.bnd_step_method = HMBlay::MethFromString("KEEP_SHAPE");
	inp2.partition = {0.0, 0.01, 0.02, 0.03, 0.04};
	inp2.bnd_step = 0.01;
	inp2.start=inp2.end=Point(0,0);
	GridGeom Ans3 = HMBlay::BuildBLayerGrid({inp2});
	save_vtk(Ans3, "t01_3.vtk");
	add_check( [&]()->bool{
		//points should lie strictly on the source contour
		auto gcont = GGeom::Info::Contour(Ans3);
		HMCont2D::Contour& inner = *gcont.roots()[0]->children[0];
		for (auto pt: inner.all_points()){
			double dist = std::get<4>(col2.coord_at(*pt));
			if (!ISZERO(dist)) {return false;}
		}
		return true;
	}(), cn);
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

	std::string cn("8-side polygon with outer layer");
	auto col1 = HMCont2D::Constructor::Circle(8, 1.3, Point(2, 2));
	inp.edges = &col1;
	GridGeom Ans1 = HMBlay::BuildBLayerGrid({inp});

	vector<double> area1;
	for(int i = 2*Ans1.n_cells()/3+1; i<Ans1.n_cells(); ++i)
		area1.push_back(Ans1.get_cell(i)->area());
	//maximum outer element is no more then two times bigger then minimum
	double mina1 = *(std::min_element(area1.begin(), area1.end()));
	double maxa1 = *(std::max_element(area1.begin(), area1.end()));
	add_check(mina1>0 && maxa1/mina1<2.1, cn);

	cn = std::string("8-side polygon with inner layer");
	inp.direction = HMBlay::DirectionFromString("INNER");
	inp.partition = {0.0, 0.02, 0.08, 0.2};
	GridGeom Ans2 = HMBlay::BuildBLayerGrid({inp});
	vector<double> area2;
	for(int i = 0; i<Ans2.n_cells()/3; ++i)
		area2.push_back(Ans2.get_cell(i)->area());
	//maximum outer element is bigger not more then two times then minimum
	double mina2 = *(std::min_element(area2.begin(), area2.end()));
	double maxa2 = *(std::max_element(area2.begin(), area2.end()));
	add_check(mina2>0 && maxa2/mina2<2.1, cn);
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
	auto col1 = HMCont2D::Constructor::ContourFromPoints(
			{0,0, 1,1, 2,1.9, 4,3.5, 7,6.1});
	inp.direction = HMBlay::DirectionFromString("LEFT");
	inp.start = Point(0,0);
	inp.end = Point(7,6.1);
	inp.edges = &col1;
	GridGeom Ans1 = HMBlay::BuildBLayerGrid({inp});
	add_check([&](){
		for (int i=0; i<Ans1.n_cells(); ++i) if (Ans1.get_cell(i)->area()<0) return false;
		for (int i=0; i<Ans1.n_points(); ++i) if (Ans1.get_point(i)->y<0.0) return false;
		return true;
	}(), cn);

	cn = std::string("5-points polyline. Right layer");
	inp.direction = HMBlay::DirectionFromString("RIGHT");
	GridGeom Ans2 = HMBlay::BuildBLayerGrid({inp});
	add_check([&](){
		for (int i=0; i<Ans2.n_cells(); ++i) if (Ans2.get_cell(i)->area()<0) return false;
		for (int i=0; i<Ans2.n_points(); ++i) if (Ans2.get_point(i)->y>6.1) return false;
		return true;
	}(), cn);

	cn = std::string("5-points polyline. Right layer. Opposite direction");
	std::swap(inp.start, inp.end);
	GridGeom Ans3 = HMBlay::BuildBLayerGrid({inp});
	add_check([&](){
		for (int i=0; i<Ans3.n_cells(); ++i) if (Ans3.get_cell(i)->area()<0) return false;
		for (int i=0; i<Ans3.n_points(); ++i) if (Ans3.get_point(i)->y<0.0) return false;
		return true;
	}(), cn);
}

void test04(){
	std::cout<<"04. Different partitions"<<std::endl;

	auto col1 = HMCont2D::Constructor::Circle(8, 10, Point(0, 0));
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
	GridGeom Ans1 = HMBlay::BuildBLayerGrid({inp1});
	add_check([&](){
		for (int i=0; i<Ans1.n_cells(); ++i) if (Ans1.get_cell(i)->area()<0) return false;
		for (int i=0; i<Ans1.n_points(); ++i) if (Ans1.get_point(i)->y<-1e-3) return false;
		bool hasp1=false, hasp2=false;
		for (int i=0; i<Ans1.n_points(); ++i){
			if (Point::dist(*Ans1.get_point(i), Point(10.4, 0))<1e-2) hasp1 = true;
			if (Point::dist(*Ans1.get_point(i),Point(-10.4, 0))<1e-2) hasp2 = true;
		}
		return true;
	}(), cn);

	cn = std::string("Two layers with same partition depth count: outer");
	HMBlay::Input inp2(inp1);
	inp2.partition = {0, 0.1, 0.2, 0.3, 0.8};
	std::swap(inp2.start, inp2.end);
	inp2.bnd_step = 0.1;
	GridGeom Ans2 = HMBlay::BuildBLayerGrid({inp1, inp2});
	add_check([&](){
		for (int i=0; i<Ans2.n_cells(); ++i) if (Ans2.get_cell(i)->area()<0) return false;
		if (Ans2.n_cells()/4 != Ans2.n_points()/5) return false;
		for (int i=0; i<Ans2.n_points(); ++i) if (Ans2.get_point(i)->y<-10-0.8-1e-12) return false;
		for (int i=0; i<Ans2.n_points(); ++i) if (Ans2.get_point(i)->y>10+0.4+1e-12) return false;
		bool hasp1=false, hasp2=false;
		for (int i=0; i<Ans2.n_points(); ++i){
			if (Point::dist(*Ans2.get_point(i), Point(10.8, 0))<1e-2) hasp1 = true;
			if (Point::dist(*Ans2.get_point(i),Point(-10.8, 0))<1e-2) hasp2 = true;
		}
		if (!hasp1 || !hasp2) return false;
		return true;
	}(), cn);

	cn = std::string("Two layers with same partition depth count: inner");
	inp1.direction = inp2.direction = HMBlay::DirectionFromString("INNER");
	GridGeom Ans3 = HMBlay::BuildBLayerGrid({inp1, inp2});
	add_check([&](){
		for (int i=0; i<Ans3.n_cells(); ++i) if (Ans3.get_cell(i)->area()<0) return false;
		if (Ans3.n_cells()/4 != Ans3.n_points()/5) return false;
		for (int i=0; i<Ans3.n_points(); ++i){
			const GridPoint* p=Ans3.get_point(i);
			if (p->y>10+1e-3 || p->y<-10-1e-3) return false;
			if (ISEQ(p->x, 0) && (p->y>-9.13 && p->y<9.45)) return false;
		}
		return true;
	}(), cn);

	cn = std::string("Two layers with different partition depth count");
	inp1.direction = inp2.direction = HMBlay::DirectionFromString("OUTER");
	inp2.partition.push_back(1.0); inp2.partition.push_back(1.2);
	GridGeom Ans4 = HMBlay::BuildBLayerGrid({inp1, inp2});
	add_check([&](){
		for (int i=0; i<Ans4.n_cells(); ++i) if (Ans4.get_cell(i)->area()<0) return false;
		int cp = 0;
		for (int i=0; i<Ans4.n_points(); ++i){
			if (fabs(Ans4.get_point(i)->y)<0.02){
				if (Ans4.get_point(i)->x >  10.81) ++cp;
				if (Ans4.get_point(i)->x < -10.81) ++cp;
			}
		}
		if (cp!=4) return false;
		return true;
	}(), cn);
}


void test05(){
	std::cout<<"05. Non-right boundary angle at the ends"<<std::endl;
	auto col1 = HMCont2D::Constructor::ContourFromPoints(
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
	GridGeom Ans1 = HMBlay::BuildBLayerGrid({inp1});
	add_check([&](){
		if (fabs(Ans1.area()-1.93)>0.01) return false;
		for (int i=0; i<Ans1.n_cells(); ++i) if (Ans1.get_cell(i)->area()<0) return false;
		return true;
	}(), cn);
	save_vtk(Ans1, "t5.vtk");
}


void test06(){
	std::cout<<"06. 1 segment sources"<<std::endl;
	auto col1 = HMCont2D::Constructor::ContourFromPoints(
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
	GridGeom Ans1 = HMBlay::BuildBLayerGrid({inp1});
	add_check([&](){
		if (Ans1.n_cells() != 1) return false;
		for (int i=0; i<Ans1.n_points(); ++i) if (Ans1.get_point(i)->y<-1e-12) return false;
		return true;
	}(), cn);

	cn = std::string("right from positive");
	inp1.direction = HMBlay::DirectionFromString("RIGHT");
	GridGeom Ans2 = HMBlay::BuildBLayerGrid({inp1});
	add_check([&](){
		if (Ans2.n_cells() != 1) return false;
		for (int i=0; i<Ans2.n_points(); ++i) if (Ans2.get_point(i)->y>1e-12) return false;
		return true;
	}(), cn);

	cn = std::string("left from negative");
	std::swap(inp1.start, inp1.end);
	inp1.direction = HMBlay::DirectionFromString("LEFT");
	GridGeom Ans3 = HMBlay::BuildBLayerGrid({inp1});
	add_check([&](){
		if (Ans3.n_cells() != 1) return false;
		for (int i=0; i<Ans3.n_points(); ++i) if (Ans3.get_point(i)->y>1e-12) return false;
		return true;
	}(), cn);

	cn = std::string("right from negative");
	inp1.direction = HMBlay::DirectionFromString("RIGHT");
	GridGeom Ans4 = HMBlay::BuildBLayerGrid({inp1});
	add_check([&](){
		if (Ans4.n_cells() != 1) return false;
		for (int i=0; i<Ans4.n_points(); ++i) if (Ans4.get_point(i)->y<-1e-12) return false;
		return true;
	}(), cn);

	cn = std::string("right from negative. fine partition");
	inp1.bnd_step_method = HMBlay::MethFromString("IGNORE_ALL");
	GridGeom Ans5 = HMBlay::BuildBLayerGrid({inp1});
	add_check([&](){
		if (Ans5.n_cells() != 20) return false;
		for (int i=0; i<Ans5.n_points(); ++i) if (Ans5.get_point(i)->y<-1e-12) return false;
		return true;
	}(), cn);
}


void test07(){
	std::cout<<"07. Corner angle algorithm basic test"<<std::endl;
	auto c1 = HMCont2D::Constructor::Circle(36, 3, Point(0, 0));
	HMCont2D::PCollection apoints;
	apoints.add_value(Point(-3, -7));
	apoints.add_value(Point(4, -7));
	HMCont2D::Contour con1 = HMCont2D::Assembler::Contour1(c1,
			HMCont2D::ECollection::FindClosestNode(c1, Point(3,0)),
			HMCont2D::ECollection::FindClosestNode(c1, Point(-3,0)));
	con1.add_value(HMCont2D::Edge(con1.last(), apoints.point(0)));
	con1.add_value(HMCont2D::Edge(apoints.point(0), apoints.point(1)));
	con1.add_value(HMCont2D::Edge(apoints.point(1), con1.first()));
	
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
	GridGeom Ans1 = HMBlay::BuildBLayerGrid({inp1});
	add_check(Ans1.n_points() == 545 && Ans1.n_cells() == 436,
			"domain with a half-cirlce");

	save_vtk(Ans1, "t7.vtk");
}

void test08(){
	std::cout<<"08 Single angle connection"<<std::endl;
	auto c1 = HMCont2D::Constructor::Circle(36, 3, Point(0, 0));
	HMCont2D::PCollection apoints;
	apoints.add_value(Point(3, -3));
	HMCont2D::Contour con1 = HMCont2D::Assembler::Contour1(c1,
			HMCont2D::ECollection::FindClosestNode(c1, Point(3,0)),
			HMCont2D::ECollection::FindClosestNode(c1, Point(0,-3)));
	con1.add_value(HMCont2D::Edge(con1.last(), apoints.point(0)));
	con1.add_value(HMCont2D::Edge(apoints.point(0), con1.first()));

	HMBlay::Input inp1;
	inp1.bnd_step_method = HMBlay::MethFromString("KEEP_SHAPE");
	inp1.direction = HMBlay::DirectionFromString("INNER");
	inp1.partition = {0, 0.02, 0.05, 0.1, 0.2};
	inp1.start = inp1.end = Point(0, 0);
	inp1.bnd_step = 0.2;
	inp1.edges = &con1;

	GridGeom Ans1 = HMBlay::BuildBLayerGrid({inp1});
	add_check(Ans1.n_points() == 565 && Ans1.n_cells() == 452, "Mesh inner");

	inp1.direction = HMBlay::DirectionFromString("OUTER");
	GridGeom Ans2 = HMBlay::BuildBLayerGrid({inp1});
	add_check(Ans2.n_points() == 575 && Ans2.n_cells() == 460, "Mesh outer");

	save_vtk(Ans1, "t8_1.vtk");
	save_vtk(Ans2, "t8_2.vtk");
}

void test09(){
	std::cout<<"09 Throw on crossing intervals"<<std::endl;
	//HMFem::Mat m;
	//for (int i=0; i<100; ++i) m.set(i,i,1);
	//for (int i=0; i<99; ++i) m.set(i, i+1, 0.3);
	//auto slv = HMFem::MatSolve::Factory(m);
}

void test10(){
	std::cout<<"10. right + re-entrant angle"<<std::endl;
	auto c = HMCont2D::Constructor::ContourFromPoints(
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
	GridGeom Ans1 = HMBlay::BuildBLayerGrid({inp1});
	add_check(Ans1.n_cells() == 232 && Ans1.n_points() == 295, "Mesh");
	save_vtk(Ans1, "t10.vtk");
};

void test11(){
	std::cout<<"11. Acute angle"<<std::endl;
	auto c = HMCont2D::Constructor::ContourFromPoints(
			{-10,0, 5,0, 0,1.5});
	HMBlay::Input inp1;
	inp1.bnd_step_method = HMBlay::MethFromString("KEEP_SHAPE");
	inp1.direction = HMBlay::DirectionFromString("INNER");
	inp1.edges = &c;
	inp1.bnd_step = 0.1;
	inp1.start = Point(0,0);
	inp1.end = Point(0,1.5);
	inp1.partition = {0, 0.01, 0.05, 0.1, 0.15, 0.2};

	GridGeom Ans1 = HMBlay::BuildBLayerGrid({inp1});
	int badskew1=0;
	for (double s: GGeom::Info::Skewness(Ans1)) {if (s>0.7) ++badskew1;}
	add_check(Ans1.n_cells() > 550 && Ans1.n_points() > 590 && badskew1==1, "Mesh1");

	inp1.partition = {0, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2};
	GridGeom Ans2 = HMBlay::BuildBLayerGrid({inp1});
	int badskew2=0;
	for (double s: GGeom::Info::Skewness(Ans2)) {if (s>0.7) ++badskew2;}
	add_check(Ans2.n_cells() > 700 && Ans2.n_points() > 720 && badskew2==1, "Mesh2");

	save_vtk(Ans1, "t11_1.vtk");
	save_vtk(Ans2, "t11_2.vtk");
};

void test12(){
	std::cout<<"12. Different kinds of reflex angles"<<std::endl;
	auto c = HMCont2D::Constructor::ContourFromPoints(
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

	GridGeom Ans1 = HMBlay::BuildBLayerGrid({inp1});
	add_check(Ans1.n_points() == 938 && Ans1.n_cells() == 785, "Mesh");
	save_vtk(Ans1, "t12.vtk");
}

void test13(){
	std::cout<<"13. Cylinder in channel meshing"<<std::endl;
	//g1
	auto g1 = GGeom::Constructor::RectGrid(Point(0,0), Point(30, 10), 90, 30);
	//g2
	HMBlay::Input inp1, inp2;
	auto c = HMCont2D::Constructor::ContourFromPoints({0, 0, 30, 0});
	inp1.bnd_step_method = HMBlay::MethFromString("KEEP_SHAPE");
	inp1.direction = HMBlay::DirectionFromString("LEFT");
	inp1.edges = &c;
	inp1.bnd_step = 0.2;
	inp1.start = Point(0,0);
	inp1.end = Point(30, 0);
	inp1.partition = {0, 0.03, 0.07, 0.12, 0.2, 0.3, 0.45, 0.7};
	GridGeom g2_1 = HMBlay::BuildBLayerGrid({inp1});
	c = HMCont2D::Constructor::ContourFromPoints({30, 10, 0, 10});
	inp1.edges = &c;
	inp1.start = Point(30,10); inp1.end = Point(0, 10);
	GridGeom g2_2 = HMBlay::BuildBLayerGrid({inp1});
	GridGeom g2 = GridGeom::sum(vector<GridGeom*> {&g2_1, &g2_2});
	//g3
	c = HMCont2D::Constructor::Circle(64, 0.5, Point(10, 5));
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
	GridGeom g3 = HMBlay::BuildBLayerGrid({inp1, inp2});
	//g4
	auto conf = HMMath::Conformal::Rect::Factory(
		{Point(10.38, 4.07), Point(15.33, 3), Point(15.33, 7), Point(10.38, 5.93)},
		{0,1,2,3});
	auto g4 = GGeom::Constructor::RectGrid(Point(0, 0), Point(conf->module(), 1), 30, 50);
	GGeom::Modify::PointModify(g4, [conf](GridPoint* p){ Point a = conf->MapToPolygon1(*p); p->set(a.x, a.y);});
	//g5
	auto g5 = GGeom::Constructor::RectGrid(Point(14.5, 2.5), Point(30, 7.5), 92, 30);

	//gcommon
	auto gr1 = GridGeom::cross_grids(&g1, &g5, 0.3, 7, false, false, CrossGridCallback::to_cout());
	auto gr2 = GridGeom::cross_grids(gr1, &g4, 0.1, 7, false, false, CrossGridCallback::to_cout());
	auto gr3 = GridGeom::cross_grids(gr2, &g2, 0.1, 7, false, false, CrossGridCallback::to_cout());
	auto gr4 = GridGeom::cross_grids(gr3, &g3, 0.2, 7, false, true, CrossGridCallback::to_cout());

	//save to file
	save_vtk(g1, "t13_g1.vtk");
	save_vtk(g2, "t13_g2.vtk");
	save_vtk(g3, "t13_g3.vtk");
	save_vtk(g4, "t13_g4.vtk");
	save_vtk(g5, "t13_g5.vtk");
	save_vtk(gr4, "t13_gres.vtk");

	delete gr1; delete gr2; delete gr3; delete gr4;
}

GridGeom* grid_minus_cont(GridGeom& g, const HMCont2D::Contour& cont){
	vector<Point> pp;
	vector<int> eds;
	for (auto p: cont.ordered_points()) pp.push_back(*p);
	pp.resize(pp.size()-1);
	for (int i=0; i<pp.size(); ++i){
		eds.push_back(i);
		eds.push_back(i+1);
	}
	eds.back() = 0;
	PointsContoursCollection col(pp, eds);
	return GridGeom::grid_minus_cont(&g, &col, false, CrossGridCallback::to_cout());
}

void test14(){
	std::cout<<"14. Four angle types"<<std::endl;

	//right angle
	auto fx = [](double t){ return 0.01*sin(2.0*M_PI*3.0*t); };
	auto fy = [](double t){ return 0.013*(-sin(2.0*M_PI*5.0*t)); };
	double h=0.01;
	double t=1-h;
	HMCont2D::PCollection pc1;
	while (t>h/2.0){
		pc1.add_value(Point(fx(t)+0.08*t, t));
		t-=h;
	}
	while (t<1+h/2.0){
		pc1.add_value(Point(t, fy(t)));
		t+=h;
	}
	pc1.add_value(Point(1,1));
	auto c1 = HMCont2D::Constructor::ContourFromPoints(pc1, true); 
	HMBlay::Input inp1;
	inp1.bnd_step_method = HMBlay::MethFromString("IGNORE_ALL");
	inp1.direction = HMBlay::DirectionFromString("LEFT");
	inp1.edges = &c1;
	inp1.bnd_step = 0.03;
	inp1.force_conformal = true;
	inp1.partition = {0, 0.01, 0.02, 0.03, 0.04, 0.053, 0.07, 0.095, 0.12};

	GridGeom g1 = HMBlay::BuildBLayerGrid({inp1});
	GridGeom basgrid1 = GGeom::Constructor::RectGrid(Point(-0.1, -0.1), Point(1.1, 1.1), 30, 30);
	GridGeom* b1 = grid_minus_cont(basgrid1, *GGeom::Info::Contour(g1).roots()[0]);
	GridGeom* gr1 = GridGeom::cross_grids(b1, &g1, 0.03, 7, false, false, CrossGridCallback::to_cout());
	GGeom::Modify::SimplifyBoundary(*gr1, M_PI/4);
	add_check(fabs(gr1->area() - 0.955398)<1e-5, "Square with curved boundary");
	save_vtk(gr1, "t14_g1.vtk");

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
	auto c2 = HMCont2D::Constructor::ContourFromPoints(pc2, true);
	inp1.edges=&c2;
	inp1.bnd_step = 0.02;
	inp1.start =inp1.end = Point(0, 0);
	GridGeom g2 = HMBlay::BuildBLayerGrid({inp1});
	GridGeom basgrid2 = GGeom::Constructor::RectGrid(Point(-2.1, -1.1), Point(2.8, 2.1), 150, 100);
	GridGeom* b2 = grid_minus_cont(basgrid2, *GGeom::Info::Contour(g2).roots()[0]);
	GridGeom* gr2 = GridGeom::cross_grids(b2, &g2, 0.03, 7, false, false, CrossGridCallback::to_cout());
	GGeom::Modify::SimplifyBoundary(*gr2, M_PI/4);
	add_check(fabs(gr2->area() - 6.26758)<1e-5, "Reentrant area with a hole");
	save_vtk(gr2, "t14_g2.vtk");

	//acute angle
	auto c3 = HMCont2D::Constructor::ContourFromPoints(
			{0,0, 5,0, 0,1.5}, true);
	HMBlay::Input inp2;
	inp2.bnd_step_method = HMBlay::MethFromString("KEEP_SHAPE");
	inp2.direction = HMBlay::DirectionFromString("INNER");
	inp2.edges = &c3;
	inp2.bnd_step = 0.07;
	inp2.force_conformal = true;
	inp2.start = inp2.end = Point(0,0);
	inp2.partition = {0, 0.02, 0.05, 0.1, 0.15, 0.2};
	GridGeom g3 = HMBlay::BuildBLayerGrid({inp2});
	GridGeom basgrid3 = GGeom::Constructor::RectGrid(Point(-0.1, -0.1), Point(4.5, 2), 80 , 40);
	GridGeom* r3 = GridGeom::cross_grids(&basgrid3, &g3, 0.03, 7, true, false, CrossGridCallback::to_cout());
	GridGeom* gr3 = grid_minus_cont(*r3, c3);
	add_check(fabs(gr3->area() - 3.75)<1e-5, "Closed tringle with acute angle");
	save_vtk(gr3, "t14_g3.vtk");
	
};

void test15(){
	std::cout<<"15. Geometry noise"<<std::endl;

	auto c1 = HMCont2D::Constructor::ContourFromPoints({0,0, 0.5, -0.03, 0.55, 0, 0.68, 0, 0.7, 0.02, 0.75, 0, 1, 0, 1.01, -0.02, 1.03, 0, 1.1, 0}); 
	HMCont2D::SaveVtk(c1, "c15.vtk");
	HMBlay::Input inp1;
	inp1.bnd_step_method = HMBlay::MethFromString("KEEP_SHAPE");
	inp1.direction = HMBlay::DirectionFromString("LEFT");
	inp1.edges = &c1;
	inp1.bnd_step = 0.05;
	inp1.start = Point(0, 0);
	inp1.end = Point(10, 0);
	inp1.partition = {0, 0.01, 0.02, 0.04, 0.08, 0.12, 0.16, 0.2};
	GridGeom g1 = HMBlay::BuildBLayerGrid({inp1});

	save_vtk(g1, "t15.vtk");
}

void test16(){
	std::cout<<"16. Impostition of a boundary grid with contour conflicts"<<std::endl;
	//1) cicle grid
	GridGeom gcirc = GGeom::Constructor::Ring(Point(0, 0), 1, 0.6, 34, 5);

	//2) get boundary from it
	HMCont2D::ContourTree gcont = GGeom::Info::Contour(gcirc);
	//3) build a boundary layer around it
	HMBlay::Input inp1;
	inp1.bnd_step_method = HMBlay::MethFromString("IGNORE_ALL");
	inp1.direction = HMBlay::DirectionFromString("LEFT");
	inp1.edges = &gcont;
	inp1.bnd_step = 0.06;
	inp1.start = Point(1, 0);
	inp1.end = Point(1, 1);
	inp1.partition = {0, 0.01, 0.02, 0.04};
	GridGeom bgrid = HMBlay::BuildBLayerGrid({inp1}); 
	add_check(fabs(GGeom::Info::Area(bgrid)-3.01008e-2)<1e-6, "IGNORE_ALL boundary partition");
	//simplify boundary tests
	GridGeom bgrid2 = GGeom::Constructor::DeepCopy(bgrid);
	GGeom::Modify::SimplifyBoundary(bgrid2, M_PI);
	add_check(fabs(GGeom::Info::Area(bgrid2)-2.82039e-2)<1e-6, "full simplification of bgrid");
	GridGeom bgrid3 = GGeom::Constructor::DeepCopy(bgrid);
	GGeom::Modify::SimplifyBoundary(bgrid3, M_PI/2.0);
	add_check(fabs(GGeom::Info::Area(bgrid3)-2.88313e-2)<1e-6, "simplification with angle=pi/2");
	GridGeom bgrid4 = GGeom::Constructor::DeepCopy(bgrid);
	GGeom::Modify::SimplifyBoundary(bgrid4, M_PI/4.0);
	add_check(fabs(GGeom::Info::Area(bgrid4)-2.98449e-2)<1e-6, "simplification with angle=pi/4");
	GridGeom bgrid5 = GGeom::Constructor::DeepCopy(bgrid);
	GGeom::Modify::SimplifyBoundary(bgrid5, 0.0);
	add_check(fabs(GGeom::Info::Area(bgrid5)-3.01008e-2)<1e-6, "simplification with angle=0");


	//impose
	GridGeom* impgrid = GridGeom::cross_grids(&gcirc, &bgrid, 0.3, 7,
			false, false, CrossGridCallback::to_cout());
	GridGeom* impgrid1 = GridGeom::cross_grids(&gcirc, &bgrid4, 0.3, 7,
			false, false, CrossGridCallback::to_cout());
	//5) see the result
	double arinit = GGeom::Info::Area(gcirc);
	add_check(ISEQ(arinit, GGeom::Info::Area(*impgrid)), "imposition with non-simplified bgrid: area");
	add_check(ISEQ(arinit, GGeom::Info::Area(*impgrid1)), "imposition with simplified bgrid: area");
	GGeom::Modify::SimplifyBoundary(*impgrid, M_PI/4.0);
	GGeom::Modify::SimplifyBoundary(*impgrid1, M_PI/4.0);
	add_check((fabs(1.99894 - GGeom::Info::Area(*impgrid))<1e-4), "area1 after simplification");
	add_check(ISEQ(arinit, GGeom::Info::Area(*impgrid1)), "area2 after simplification. Keep all cell non-degenerate");
	
	auto skew1 = GGeom::Info::Skewness(*impgrid);
	auto skew2 = GGeom::Info::Skewness(*impgrid1);
	double maxskew1 = *std::max_element(skew1.begin(), skew1.end());
	double maxskew2 = *std::max_element(skew2.begin(), skew2.end());
	vector<double> badskew1, badskew2;
	std::copy_if(skew1.begin(), skew1.end(), std::back_inserter(badskew1), [](double a){return a>0.8;});
	std::copy_if(skew2.begin(), skew2.end(), std::back_inserter(badskew2), [](double a){return a>0.8;});
	add_check(badskew1.size() == 0 && badskew2.size() == 4, "skewness check");

	save_vtk(impgrid, "t16_1.vtk");
	save_vtk(impgrid1, "t16_2.vtk");
	delete impgrid;
	delete impgrid1;
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

	HMCont2D::Container<HMCont2D::Contour>
	cont = HMCont2D::Constructor::ContourFromPoints(inppts, false);

	HMBlay::Input inp1;
	inp1.bnd_step_method = HMBlay::MethFromString("NO");
	inp1.direction = HMBlay::DirectionFromString("LEFT");
	inp1.edges = &cont;
	inp1.start = Point(0, 0);
	inp1.end = Point(2, 0);
	inp1.partition = {0, 0.01, 0.02, 0.04, 0.08, 0.14, 0.2};
	GridGeom bgrid = HMBlay::BuildBLayerGrid({inp1}); 
	HMCont2D::SaveVtk(cont, "t17_1.vtk");
	save_vtk(bgrid, "t17_2.vtk");
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
	

	if (FAILED_CHECKS ==1){
		std::cout<<FAILED_CHECKS<<" test failed <<<<<<<<<<<<<<<<<<<"<<std::endl;
	} else if (FAILED_CHECKS > 1) {
		std::cout<<FAILED_CHECKS<<" tests failed <<<<<<<<<<<<<<<<<<<"<<std::endl;
	} else {
		std::cout<<"All tests passed"<<std::endl;
	}
	std::cout<<"DONE"<<std::endl;
}
