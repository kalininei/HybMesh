#ifndef HMCONT2D_CONTCLIPPING_HPP
#define HMCONT2D_CONTCLIPPING_HPP

#include "primitives2d.hpp"
#include "contour.hpp"
#include "contour_tree.hpp"

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

}

namespace Algos{

//Offset contour
enum class OffsetTp{
	//rounding at corners
	RC_CLOSED_POLY,
	RC_OPEN_ROUND,
	RC_OPEN_BUTT,
	//square corners
	SC_CLOSED_POLY,
};

//takes into account direction of source and sign of delta:
//offsets to the right side if delta>0.
Tree Offset(const EdgeData& source, double delta, OffsetTp tp);
//forces singly connected output contour. tp = CLOSED_POLY or OPEN_ROUND
EdgeData Offset1(const EdgeData& source, double delta);

}

}}

#endif
