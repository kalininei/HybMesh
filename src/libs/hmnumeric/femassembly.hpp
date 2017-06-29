#ifndef HYBMESH_HMFEM_FEMASSEMBLY_HPP
#define HYBMESH_HMFEM_FEMASSEMBLY_HPP

#include "spmat.hpp"
#include "femgrid43.hpp"

namespace HMFem{

namespace Assemble{

//== grad(p_i) . grad(p_j)
shared_ptr<HMMath::Mat> PureLaplace(const HM2D::GridData& grid);

//== p_i * p_j
shared_ptr<HMMath::Mat> FullMass(const HM2D::GridData& grid);

//== p_i
vector<double> LumpMass(const HM2D::GridData& grid);

//== dp_j /dx * p_i 
shared_ptr<HMMath::Mat> DDx(const HM2D::GridData& grid);

//== dp_j /dy * p_i 
shared_ptr<HMMath::Mat> DDy(const HM2D::GridData& grid);

}//Assemble










}




#endif
