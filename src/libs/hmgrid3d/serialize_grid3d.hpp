#ifndef SERIALIZE_GRID3D_HPP
#define SERIALIZE_GRID3D_HPP
#include "primitives_grid3d.hpp"
#include <type_traits>
#include "hmcallback.hpp"

namespace HMGrid3D{

//only linear data is saved here
struct SimpleSerialize{
	//primitives sizes
	int n_vert, n_faces, n_edges, n_cells;

	//data
	std::vector<double> vert;  //x0, y0, z0, x1, y1, z1, ...
	std::vector<int> edges;    //edge0_start, edge0_end, edge1_start, edge1_end, ...
	std::vector<int> faces;    //face0_num_edges, face0_edge0, ..., face0_leftcell, face0_rightcell, ...
	std::vector<int> cells;    //cell0_num_faces, cell0_face0, cell0_face1, ....
	std::vector<int> bnd;      //btype0 num_faces_in_btype0 btype0_face0, ....

	//access to cells, faces in 'cells', 'faces' vectors
	vector<int> icells, ifaces;

	//fills everything from: vert, edges, faces, bnd.
	void supplement();

	//====== additional connectivity tables
	vector<int> face_vertex(int num_face) const;
	vector<vector<int>> face_vertex() const;
};


//grid with serialized data
struct SGrid: public SimpleSerialize, public GridData{

	//====== constructing procedures
	//fills SimpleSerialize on the basis of GridData
	void actualize_serial_data();
	//fills GridData vectors from SimpleSerialize Data
	void actualize_data();

	//====== methods
	void set_btype(std::function<int(Vertex, int)> func);
	void renumber_by_cells();
};


}

#endif
