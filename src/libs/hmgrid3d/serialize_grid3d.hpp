#ifndef SERIALIZE_GRID3D_HPP
#define SERIALIZE_GRID3D_HPP
#include "primitives_grid3d.hpp"
#include "surface_grid3d.hpp"
#include <type_traits>
#include "hmcallback.hpp"

namespace HMGrid3D{

struct SimpleSerializeSurface{
	SimpleSerializeSurface():cache(this){}
	void supplement();

	size_t n_vert, n_edges, n_faces;
	std::vector<double> vert;  //x0, y0, z0, x1, y1, z1, ...
	std::vector<int> edges;    //edge0_start, edge0_end, edge1_start, edge1_end, ...
	std::vector<std::vector<int>> face_edge;     //face-edge connectivity
	std::vector<int> btypes;                     //boundary type for each face

	//====== additional connectivity tables
	const vector<int>& face_vertex(int num_face) const { return cache.face_vertex()[num_face]; }
	const vector<vector<int>>& face_vertex() const { return cache.face_vertex(); }

	//cached data
	struct Cache{
		Cache(const SimpleSerializeSurface* parent): parent(parent){};
		const SimpleSerializeSurface* parent;

		shared_ptr<vector<vector<int>>> _face_vertex;
		vector<vector<int>>& face_vertex();

		void clear(){
			_face_vertex.reset();
		};
	};
	mutable Cache cache;
};

//surface with serialized data
struct SSurface: public SimpleSerializeSurface, public Surface{
	SSurface(){};
	SSurface(SimpleSerializeSurface&& ss) noexcept;
	SSurface(HMGrid3D::Surface&& srf) noexcept;

	//====== constructing procedures
	//fills SimpleSerialize on the basis of Surface
	void actualize_serial_data() noexcept;
	//fills Surface vectors from SimpleSerialize Data
	void actualize_data() noexcept;
};

struct SimpleSerialize{
	SimpleSerialize(): cache(this){}
	//fills everything from: vert, edges, face_edge, face_cell, btypes
	void supplement();

	//primitives sizes
	size_t n_vert, n_faces, n_edges, n_cells;

	//data
	std::vector<double> vert;                 //x0, y0, z0, x1, y1, z1, ...
	std::vector<int> edges;                   //edge0_start, edge0_end, edge1_start, edge1_end, ...
	std::vector<std::vector<int>> face_edge;  //face-edge connectivity
	std::vector<int> face_cell;               //face-cell connectivity
	std::vector<int> btypes;                  //boundary type for each face

	//====== extracts surface
	SimpleSerializeSurface serialized_surface() const;

	//====== additional connectivity tables
	const vector<int>& face_vertex(int num_face) const { return cache.face_vertex()[num_face]; }
	const vector<vector<int>>& face_vertex() const { return cache.face_vertex(); }
	const vector<vector<int>>& cell_face() const { return cache.cell_face(); }
	const vector<vector<int>>& cell_vertex() const { return cache.cell_vertex(); }
	const vector<int>& bvert() const { return cache.bvert(); }
	const vector<int>& bedges() const { return cache.bedges(); }
	const vector<int>& bfaces() const { return cache.bfaces(); }

	//cached data
	struct Cache{
		Cache(const SimpleSerialize* parent): parent(parent){};
		const SimpleSerialize* parent;

		shared_ptr<vector<vector<int>>> _face_vertex;
		vector<vector<int>>& face_vertex();

		shared_ptr<vector<vector<int>>> _cell_face;
		vector<vector<int>>& cell_face();

		shared_ptr<vector<vector<int>>> _cell_vertex;
		vector<vector<int>>& cell_vertex();

		shared_ptr<vector<int>> _bfaces;
		vector<int>& bfaces();

		shared_ptr<vector<int>> _bedges;
		vector<int>& bedges();

		shared_ptr<vector<int>> _bvert;
		vector<int>& bvert();

		void clear(){
			_face_vertex.reset();
			_bfaces.reset();
			_bedges.reset();
			_bvert.reset();
		};
	};
	mutable Cache cache;
};


//grid with serialized data
struct SGrid: public SimpleSerialize, public GridData{
	SGrid(){}
	SGrid(GridData&& gd) noexcept: GridData(std::move(gd)){
		actualize_serial_data();
	}
	explicit SGrid(const GridData& gd): GridData(gd){
		actualize_serial_data();
	}

	//====== constructing procedures
	//fills SimpleSerialize on the basis of GridData
	void actualize_serial_data() noexcept;
	//fills GridData vectors from SimpleSerialize Data
	void actualize_data() noexcept;

	//====== methods
	void set_btype(std::function<int(Vertex, int)> func);
	void renumber_by_cells();
};


}

#endif
