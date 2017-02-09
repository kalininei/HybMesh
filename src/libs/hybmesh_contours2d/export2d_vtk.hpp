#ifndef HYBMESH_VTK_EXPORT2D
#define HYBMESH_VTK_EXPORT2D

#include "primitives2d.hpp"

namespace HM2D{ namespace Export{

void ContourVTK(const EdgeData& data, std::string fn);

void VerticesVTK(const VertexData& data, std::string fn);


void GridVTK(const GridData& g, std::string fn);

void BoundaryVTK(const GridData& g, std::string fn);

}}
#endif
