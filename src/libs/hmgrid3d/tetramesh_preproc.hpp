#ifndef HMGRID3D_TETRAMESH_PREPROCESSOR_HPP
#define HMGRID3D_TETRAMESH_PREPROCESSOR_HPP
#include "surface_tree.hpp"

namespace HM3D{namespace Mesher{

//splits surfaces by relatively smooth subareas,
//builds pyramids on all non-triangle surface cells.
//input tree should be correctly directed and have nesting level <=1
struct SurfacePreprocess{
	SurfacePreprocess(const HM3D::Surface::Tree& tree, double split_angle);
	void Restore(GridData& g, const VertexData& gvert, const VertexData& svert);

	//======== data filled by constructor
	//data for splitted triangulated surface
	vector<vector<FaceData>> decomposed_surfs;
	EdgeData ae;
	VertexData av;
	vector<vector<vector<EdgeData>>> decomposed_edges;
	
	//boundary grid
	GridData bnd_grid;
private:
	std::map<FaceData*, Vect3> surfs_rnormals;

	//constructor subroutings
	void supplement_from_decomposed_surfs();
	void assemble_bnd_grid();
	//data for splitted surface before triangulation
	vector<vector<FaceData>> decomposed_surfs1;
	FaceData allfaces1;

	void merge_with_bnd(GridData& tar, const vector<int>& tar_points, const vector<int>& bnd_points);
};



}}



#endif
