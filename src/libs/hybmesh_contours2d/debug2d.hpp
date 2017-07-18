#ifndef HYBMESH_CONTOURS2D_DEBUG_HPP
#define HYBMESH_CONTOURS2D_DEBUG_HPP

#ifndef NDEBUG

#include "primitives2d.hpp"
#include "contour.hpp"
#include "contour_tree.hpp"
#include "hmdebug.hpp"
#include "hmtimer.hpp"

namespace HM2D{

struct Debug: public HMDebug{
	static void info_contour(const EdgeData& c);
	static void info_tree(const Contour::Tree& c);
	static void info_tree_node(const Contour::Tree::TNode& c);
	static void info_cell(const GridData& g, int icell);
	static void numer_all(const GridData& g);

	//contour in geogebra script format
	static void geogebra_contour(const EdgeData& c);
	static void geogebra_ecollection(const EdgeData& ecol);
	static void geogebra_tree(const Contour::Tree& c);
	static void geogebra_box(const BoundingBox& c);

	//to vtk format
	static void save_grid_vtk(const GridData& c);
	static void save_grid_cell_data_vtk(const GridData& c, const vector<double>& dt);
	static void save_grid_vertex_data_vtk(const GridData& c, const vector<double>& dt);
	static void save_cells_vtk(const CellData& c);
	static void save_cells_vtk_id(const CellData& c);

	static void save_edges_vtk(const EdgeData& c);
	static void save_edges_vtk(const Contour::Tree& c);
	static void save_edges_vtk(const EdgeData& c, const std::map<Point*, double>& data);
	static void save_edges_vtk_id(const EdgeData& c);

	static void save_vertices_vtk(const VertexData& c);
	static void save_vertices_vtk(const vector<Point>& c);
	static void save_vertices_vtk_id(const VertexData& c);
	static void save_vertices_vtk(const EdgeData& c);
	static void save_vertices_vtk_id(const EdgeData& c);

	//extract primitives without id change
	static VertexData allvertices(const EdgeData& c);
	static EdgeData alledges(const CellData& c);


	static double hash(const VertexData& g){
		double ret = g.size();
		int k = 0;
		for (auto& v: g)  ret += 2.4 * ( ++k % 231) * v->x;
		for (auto& v: g)  ret -= 3.1 * ( ++k % 123) * v->y;
		return ret;
	}
	static double hash(const EdgeData& g){
		double ret = g.size();
		int k = 0;
		for (auto& v: g){
			ret += (++k % 2) * (v->vertices[0]->x * 0.5 - v->vertices[1]->y * 0.4);
			ret += (++k % 3) * (v->vertices[1]->x * 0.2 - v->vertices[0]->y * 0.7);
		}
		return ret;
	}
	static double hash(const CellData& g){
		double ret = g.size();
		int k =0;
		for (auto& c: g){
			ret += (++k % 10 - 5) * hash(c->edges);
		}
		return ret;
	}
	static double hash(const GridData& g){
		return hash(g.vcells) - hash(g.vedges) + 2. * hash(g.vvert);
	}

};

}

#endif
#endif
