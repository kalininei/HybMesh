#ifndef HMCONT2D_CONTCLIPPING_HPP
#define HMCONT2D_CONTCLIPPING_HPP

#include "primitives2d.hpp"
#include "contour.hpp"
#include "tree.hpp"

namespace HM2D{ namespace Contour{ namespace Clip{

typedef HM2D::Contour::Tree TRet;
typedef HM2D::Contour::Tree ETree;
typedef HM2D::EdgeData ECont;

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
//!!! does simplifications keeping corner and boundary type transition points
void Heal(TRet& c1);

}}}

#endif
