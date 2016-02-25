#ifndef HMCONT2D_CONTCLIPPING_HPP
#define HMCONT2D_CONTCLIPPING_HPP

#include "hybmesh_contours2d.hpp"

namespace HMCont2D{ namespace Clip{

typedef HMCont2D::Container<HMCont2D::ContourTree> TRet;
typedef HMCont2D::Container<HMCont2D::ExtendedTree> TExRet;
typedef HMCont2D::Contour ECont;
typedef HMCont2D::ContourTree ETree;

//two contours. direction is not taken into account
TRet Intersection(const ECont& c1, const ECont& c2);
TRet Union(const ECont& c1, const ECont& c2);
TRet Difference(const ECont& c1, const ECont& c2);

//tree and contour: contours direction is not taken into account
TRet Intersection(const ETree& c1, const ECont& c2);
TRet Union(const ETree& c1, const ECont& c2);
TRet Difference(const ETree& c1, const ECont& c2);
TRet Difference(const ECont& c1, const ETree& c2);

//two trees
TRet Intersection(const ETree& c1, const ETree& c2);
TRet Union(const ETree& c1, const ETree& c2);
TRet Difference(const ETree& c1, const ETree& c2);
TRet XOR(const ETree& c1, const ETree& c2);

//multiple contours operation.
//Direction of each contour is not taken into account.
TRet Intersection(const vector<ECont>& cont);
TRet Union(const vector<ECont>& cont);
TRet Difference(const ETree& c1, const vector<ECont>& cont);

//procedure to remove all bad entries which 
//may appear as a result of boolean operations
void Heal(TRet& c1);

}}

#endif
