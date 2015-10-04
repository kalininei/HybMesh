#ifndef HYBMESH_HMFEM_FEMASSEMBLY_HPP
#define HYBMESH_HMFEM_FEMASSEMBLY_HPP

#include "femmat.hpp"

namespace HMFem{

namespace Assemble{

//== grad(p_i) . grad(p_j)
shared_ptr<Mat> PureLaplas(const Grid43& grid);

//== p_i * p_j
shared_ptr<Mat> FullMass(const Grid43& grid);

//== p_i
vector<double> LumpMass(const Grid43& grid);

//== dp_j /dx * p_i 
shared_ptr<Mat> DDx(const Grid43& grid);

//== dp_j /dy * p_i 
shared_ptr<Mat> DDy(const Grid43& grid);



}//Assemble










}




#endif
