#ifndef HMGRID3D_PYRAMID_LAYER_HPP
#define HMGRID3D_PYRAMID_LAYER_HPP

#include "primitives3d.hpp"

namespace HM3D{

//Faces should have matched direction.
//Pyramids will be built to the left side of faces
//If non3only is true then triangle faces will be ignored until they form acute angle
//	with neighbouring cell.
//Resulting grid cells number will equal number of input faces.
//Input faces will be written as first faces of resulting cells.
//For non3only=true triangle faces without pyramid will form
//	single face cell with blank face-to-cell connectivity
//Inner boundary pyramid faces will have no left adjacent cell.
GridData BuildPyramidLayer(
	const FaceData& faces,
	bool non3only=true,
	double merge_angle=60);


}
#endif
