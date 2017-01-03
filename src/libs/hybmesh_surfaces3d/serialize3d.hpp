#ifndef SERIALIZE_GRID3D_HPP
#define SERIALIZE_GRID3D_HPP
#include "primitives3d.hpp"

namespace HM3D{ namespace Ser {
class Surface{
	struct Cache;
	mutable std::unique_ptr<Cache> cache;
	void empty_cache() const;
public:
	FaceData surface;
	//constructor
	Surface();
	explicit Surface(HM3D::FaceData&& srf) noexcept;
	explicit Surface(const HM3D::FaceData& srf);
	~Surface();

	//this should be called after all this->grid changes
	void reset_geometry() const { empty_cache(); }

	int n_vert() const;
	int n_edges() const;
	int n_faces() const { return surface.size(); }

	//======== main connectivity
	const vector<double>& vert() const;            //x0, y0, z0, x1, y1, z1, ...
	const vector<int>& edge_vert() const;          //edge0_start, edge0_end, edge1_start, edge1_end, ...
	const vector<vector<int>>& face_edge() const;  //face-edge connectivity
	const vector<int>& btypes() const;             //boundary type for each face

	//====== additional connectivity tables
	const vector<int>& face_vertex(int num_face) const;
	const vector<vector<int>>& face_vertex() const;

	//======= methods
	void fill_from_serial(const vector<double>& vert,
			const vector<int>& edgevert,
			const vector<vector<int>>& faceedge,
			const vector<int>& btypes);
};

class Grid{
	struct Cache;
	mutable std::unique_ptr<Cache> cache;
	void empty_cache() const;
public:
	HM3D::GridData grid;

	//constructor
	Grid();
	explicit Grid(HM3D::GridData&& grid) noexcept;
	explicit Grid(const HM3D::GridData& grid);
	~Grid();

	//this should be called after all this->grid changes
	void reset_geometry() const { empty_cache(); }

	int n_vert() const { return grid.vvert.size(); }
	int n_edges() const { return grid.vedges.size(); }
	int n_faces() const { return grid.vfaces.size(); }
	int n_cells() const { return grid.vcells.size(); }

	//======== main connectivity
	const vector<double>& vert() const;           //x0, y0, z0, x1, y1, z1, ...
	const vector<int>& edge_vert() const;         //edge0_pstart, edge0_pend, edge1_pstart, edge1_pend, ..
	const vector<vector<int>>& face_edge() const; //face->edge connectivity
	const vector<int>& face_cell() const;         //face_leftcell, face_rightcell
	const vector<int>& btypes() const;            //boundary type for each face

	//====== additional connectivity tables
	const vector<int>& face_vertex(int num_face) const;
	const vector<vector<int>>& face_vertex() const;
	const vector<vector<int>>& cell_face() const;
	const vector<vector<int>>& cell_vertex() const;
	const vector<int>& bvert() const;
	const vector<int>& bedges() const;
	const vector<int>& bfaces() const;

	//====== methods
	void set_btype(std::function<int(Vertex, int)> func);
	void renumber_by_cells();
	void fill_from_serial(const vector<double>& vert,
			const vector<int>& edgevert,
			const vector<vector<int>>& faceedge,
			const vector<int>& facecell,
			const vector<int>& btypes);
};

}}

#endif
