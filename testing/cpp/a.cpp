#include "../../build/bindings/cpp/Hybmesh.hpp"
#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>

typedef Hybmesh::Point2 P2;
typedef Hybmesh::Point3 P3;

using namespace std;

void check_cond(bool cond){
	if (!cond) throw std::runtime_error("check failed");
}
void check_dims(const vector<int>& v1, const vector<int>& v2){
	if (v1.size() != v2.size()) throw std::runtime_error("check failed");
	for (size_t i=0; i<v2.size(); ++i){
		if (v1[i] != v2[i]) throw std::runtime_error("check failed");
	}
}

int callback_good(const std::string& s1, const std::string& s2, double p1, double p2){
	std::cout<<s1<<" "<<p1<<" --- "<<s2<<" "<<p2<<std::endl;
	return 0;
}

int callback_bad(const std::string& s1, const std::string& s2, double p1, double p2){
	std::cout<<s1<<" "<<p1<<" --- "<<s2<<" "<<p2<<std::endl;
	return p1 > 0.7 ? 1 : 0;
}

void pings(){
	Hybmesh hm("../../src/py/hybmesh.py");

	auto g1 = hm.add_unf_rect_grid1({1, 2, 3}, {2, 3, 4});
	vector<int> dims = g1.dims();
	check_dims(dims, {9, 12, 4});
	auto c1 = hm.grid_bnd_to_contour(g1, true);
	auto c2 = hm.grid_bnd_to_contour(g1, false);
	check_dims(c1.dims(), {4, 4});
	check_dims(c2.dims(), {8, 8});
	vector<P2> pv(4);
	pv[0] = P2(0, 0);
	pv[1] = P2(1, 2);
	pv[2] = P2(3, 1);
	pv[3] = P2(2, 0);
	auto c3 = hm.create_contour(pv, {0, 0, 1});
	auto c4 = hm.create_contour(pv, {3, 2, 1});
	check_dims(c3.dims(), {4, 3});
	check_dims(c4.dims(), {4, 3});
	auto c5 = hm.create_spline_contour(pv, {0, 0, 1});
	auto c6 = hm.create_spline_contour(pv, {3, 2, 1});
	check_dims(c5.dims(), {101, 100});
	check_dims(c6.dims(), {101, 100});
	auto c7 = hm.extract_subcontours(c6, pv);
	check_cond(c7.size() == 3);
	auto c8 = hm.add_rect_contour({1, 1}, {10, 10}, {0, 1, 2, 3});
	check_dims(c8.dims(), {4, 4});
	auto c9 = hm.add_circ_contour({5., 5.}, 10, 32);
	check_dims(c9.dims(), {32, 32});
	double ar9 = c9.domain_area();
	check_cond(fabs(ar9 - M_PI*100) < 5);
	auto c10 = hm.simplify_contour(c2);
	check_dims(c10.dims(), {4, 4});
	auto c11 = hm.simplify_contour(c9, 91);
	std::cout<<c11.dims()[0]<<" "<<c11.dims()[1]<<std::endl;
	check_dims(c11.dims(), {4, 4});
	auto c12 = hm.unite_contours({c8, c9});
	check_dims(c12.dims(), {36, 36});
	auto c13 = hm.decompose_contour(c12);
	check_cond(c13.size() == 2);
	auto c14 = hm.clip_domain(c13[0], c13[1], "difference");
	auto c15 = hm.add_rect_contour({0, 0}, {1, 1});
	auto c16 = hm.add_rect_contour({10, 10}, {11, 11});
	auto c17 = hm.clip_domain(c15, c16, "intersection");
	auto c18 = hm.clip_domain(c13[1], c13[1], "difference");
	check_dims(c17.dims(), {0, 0});
	check_dims(c18.dims(), {0, 0});
	auto c19 = hm.partition_contour_const(c15, 0.01);
	check_dims(c19.dims(), {400, 400});
	auto c20 = hm.partition_contour_const(c15, 0.01);
	check_dims(c19.dims(), {400, 400});
	auto c21 = hm.partition_contour_ref_points(c15, {0.01, 0.5},
		{P2(0, 0), P2(1, 1)},
		30, true, -1, {}, {}, P2(0, 0), P2(1, 1));
	check_dims(c21.dims(), {18, 18});
	auto c22 = hm.partition_contour_ref_lengths(c15, {0.01, 0.5},
		{0, 2},
		30, true, -1, {}, {}, P2(0, 0), P2(1, 1));
	check_dims(c21.dims(), c22.dims());
	auto c23 = hm.matched_partition(c15, 0.1, 0.5, {},
			{0.01}, {P2(0.1, 0.1)});
	check_dims(c23.dims(), {74, 74});
	auto c24 = hm.partition_segment(0.5, 1.5, 0.1, 0.5, {0.01, 1.0});
	check_cond(c24.size() == 5 && c24[0] == 0.5);

	auto c25 = hm.create_contour({P2(0, 0), P2(1, 0), P2(2, 0)});
	auto c26 = hm.create_contour({P2(2, 0.1), P2(3, 0.1), P2(4, 1)});
	auto c27 = hm.connect_subcontours({c25, c26}, {});
	auto c28 = hm.connect_subcontours({c25, c26}, {}, "yes");
	auto c29 = hm.connect_subcontours({c25, c26}, {}, "force");
	check_dims(c27.dims(), {5, 4});
	check_dims(c28.dims(), {4, 4});
	check_dims(c29.dims(), {5, 5});

	auto c30 = hm.add_unf_rect_grid1({1, 2, 3}, {4, 5, 6, 7});
	check_dims(c30.dims(), {12, 17, 6});
	auto c31 = hm.add_unf_circ_grid(P2(0, 0), 10, 32, 5, 1.2, false);
	check_dims(c31.dims(), {160, 288, 129});
	auto c32 = hm.add_unf_ring_grid(P2(1, 1), 3, 7, 16, 2);
	check_dims(c32.dims(), {48, 80, 32});
	auto c33 = hm.add_unf_hex_grid_in_hex(P2(-1, -1), 8, 1);
	auto c34 = hm.add_unf_hex_grid_in_rect(P2(0, 0), P2(2, 3), 1);
	auto c35 = hm.add_triangle_grid(P2(0, 0), P2(1, 0), P2(0, 1), 3, {1, 2, 3});
	check_dims(c33.dims(), {216, 306, 91});
	check_dims(c34.dims(), {28, 35, 8}); 
	check_dims(c35.dims(), {10, 15, 6}); 

	auto c36 = hm.create_contour({P2(0, 0), P2(0, 1)});
	auto c37 = hm.partition_contour_const(c36, 0.1);
	auto c38 = hm.create_contour({P2(0, 0), P2(2, -0.1)});
	auto c39 = hm.partition_contour_const(c38, 0.2);
	auto c40 = hm.add_custom_rect_grid("linear", c37, c39);
	auto c41 = hm.add_custom_rect_grid_htfi(c37, c39,
			Hybmesh::Contour2D::None(),
			Hybmesh::Contour2D::None(),
			{1, 1, 1, 0.8});
	hm.assign_callback(callback_bad);
	try{
		auto c42 = hm.add_custom_rect_grid("orthogonal", c37, c39);
	} catch (Hybmesh::EUserInterrupt& e){
		std::cout<<"Interrupt catched"<<std::endl;
	}
	hm.assign_callback(callback_good);
	auto c42 = hm.add_custom_rect_grid("orthogonal", c37, c39);
	hm.reset_callback();
	check_dims(c40.dims(), {121, 220, 100});
	check_dims(c40.dims(), c41.dims());
	check_dims(c40.dims(), c42.dims());
	try{
		auto c43_ = hm.add_circ_rect_grid(P2(0, 0), -1, 0.05);
	} catch (Hybmesh::ERuntimeError& e){
		std::cout<<"Runtime error catched:"<<std::endl;
		std::cout<<e.what()<<std::endl;
	}
	auto c43 = hm.add_circ_rect_grid(P2(0, 0), 1, 0.05);
	auto c44 = hm.add_circ_rect_grid(P2(0, 0), 1, 0.05, 1., 1., "orthogonal_rect");
	check_dims(c43.dims(), c44.dims());
	auto c45 = hm.stripe(c39, {0, 0.01, 0.02, 0.05});
	auto c46 = hm.triangulate_domain(c45, {}, {}, {}, "3");
	auto c47 = hm.grid_bnd_to_contour(c45);
	auto c48 = hm.grid_bnd_to_contour(c46);
	check_cond(fabs(c47.domain_area() - c48.domain_area())<1e-8);
	auto c49 = hm.add_rect_contour(P2(0, 0), P2(1, 1));
	auto c50 = hm.partition_contour_const(c49, 0.1);
	auto c51 = hm.pebi_fill(c50, {}, {0.5}, {P2(0.5, 0.5)});
	auto c52 = hm.simple_boundary_grid(c51, {0, 0.01, 0.02});
	auto c53 = hm.simple_boundary_grid(c51, {0, 0.01, 0.02}, "right");
	check_dims(c52.dims(), c53.dims());
	auto c54 = hm.exclude_contours(c53, {c51}, "inner");
	auto c55 = hm.exclude_contours(c53, {c51}, "outer");
	check_dims(c54.dims(), c53.dims());
	check_dims(c55.dims(), {0, 0, 0});
	std::vector<double> c56x, c57x;
	for (int i=0; i<11; ++i){
		c56x.push_back(i/10.);
		c57x.push_back(0.3 + (double)i/30.);
	}
	auto c56 = hm.add_unf_rect_grid1(c56x, c56x);
	auto c57 = hm.add_unf_rect_grid1(c57x, c57x);
	auto c58 = hm.unite_grids1(c56, c57, 0.1);
	hm.stdout_verbosity(3);
	auto c59 = hm.map_grid(c58, c57, {P2(0, 0)}, {P2(0.3, 0.3)});
	hm.stdout_verbosity(0);
	hm.heal_grid(c59, 30, 30);
	auto c60 = hm.exclude_contours(c58, {c57}, "inner");
	auto c61 = hm.extrude_grid(c60, {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3});
	auto c62 = hm.grid3_bnd_to_surface(c61);
	hm.assign_callback(callback_good);
	auto c63 = hm.tetrahedral_fill({c62});
	hm.reset_callback();
	auto c64 = hm.revolve_grid(c60, P2(0, 0), P2(0, 1),
			{0, 45, 90, 180});
	auto c65 = hm.revolve_grid(c60, P2(0, 0), P2(0, 1),
			{180, 270, 360});
	hm.export3d_grid_vtk({c64}, "c64.vtk");
	hm.export3d_grid_vtk({c65}, "c65.vtk");
	auto c66 = hm.merge_grids3(c64, c65);
	hm.export3d_grid_vtk({c66}, "c66.vtk");
	auto c67 = c66.deepcopy();
	c67.scale(20, 20, 20, P3(0, 0, 0));
	auto c68 = hm.grid3_bnd_to_surface(c66);
	auto c69 = hm.grid3_bnd_to_surface(c67);
	check_cond(fabs(c68.domain_volume()-125*c69.domain_volume())<1e-8);
	c60.rotate(20, P2(0, 0));
	c69.free();

	hm.add_boundary_type(1, "bnd1");
	hm.add_boundary_type(2, "bnd2");
	auto c70 = hm.add_unf_rect_grid(P2(0, 0), P2(1, 1), 10, 10, {0, 1, 2, 3});
	hm.export_grid_msh(c70, "out.msh", {0}, {2}, {false});

	auto c71 = c70.skewness(-1);
	auto c72 = c70.skewness(0.1);
	check_cond(c71.size() == 100 && c72.size() == 0);

	c70.set_btypes_all(22);
	hm.export_contour_vtk(c70, "out.vtk");
	auto c73 = c70.raw_vertices();
	check_cond(c73.size() == 242);
	auto c74 = c70.raw_tab("bnd_bt");
	check_cond(c74.size() == 80);
	auto c75 = c70.raw_tab("bnd");
	check_cond(c75.size() == 40);
	auto c76 = c70.raw_tab("cell_edge");
	check_cond(c76.size() == 400);
	auto c77 = c70.raw_tab("edge_cell");
	check_cond(c77.size() == 440);
	auto c78 = c70.raw_tab("bt");
	check_cond(c78.size() == 220 && std::all_of(c78.begin(), c78.end(),
	                                           [](int i){ return i==22 || i==0;}));
	auto c79 = c70.raw_tab("cell_dim");
	check_cond(c79.size() == 100 && std::all_of(c79.begin(), c79.end(),
	                                           [](int i){ return i==4;}));
	auto c80 = c70.raw_tab("cell_vert");
	check_cond(c80.size() == 400);
	c70.set_btypes(1, std::vector<int>(c75.begin(), c75.begin()+3));
	hm.export_grid_tecplot(c70, "out.dat");

	auto c81 = hm.extrude_grid(c70, {0, 1, 2});
	auto c82 = hm.grid3_bnd_to_surface(c81);
	auto c83 = c81.raw_vertices();
	auto c84 = c82.raw_vertices();
	check_cond(c83.size() == 121*3*3);
	check_cond(c84.size() == (121*2 + 40)*3);
	auto c85 = c81.raw_tab("cell_vert");
	check_cond(c85.size() == 200*8);
	auto c86 = c82.raw_tab("face_vert");
	check_cond(c86.size() == (200+80)*4);

	auto c87 = hm.grid_bnd_to_contour(c70);
	auto c88 = c87.raw_vertices();
	check_cond(c88.size() == 12);
	auto c89 = c87.raw_tab("bt");
	check_cond(c89.size() == 6);

	auto c90 = hm.add_rect_contour(P2(1, 1.2), P2(2, 2), {0, 1, 2, 3});
	auto c91 = hm.add_rect_contour(P2(0, 0), P2(1.5, 1.5), {0, 1, 2, 3});
	auto c92 = hm.unite_contours({c90, c91});
	auto c93 = hm.decompose_contour(c92);
	check_cond(c93.size() == 3);
	auto c94 = hm.pick_contour(P2(-1, -1), c93);
	c94.get_point(P2(-1, -1));
	c94.get_point(P2::None(), P2(0.5, -1));
	c94.get_point(P2::None(), P2::None(), P2(1.0, -1));

	hm.export_all_hmd("out.hmd", "ascii");
	hm.remove_all();
	hm.import_grid_hmg("out.hmd", "Grid2D_1");
	hm.import3d_grid_hmg("out.hmd", "Grid3D_1");
	hm.remove_all();
}

int main(){
	//pinging all functions
	pings();
	std::cout<<"Done"<<std::endl;
}
