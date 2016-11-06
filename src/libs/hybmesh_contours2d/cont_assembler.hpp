#ifndef HMCONT2D_CONT_ASSEMBLER_HPP
#define HMCONT2D_CONT_ASSEMBLER_HPP

#include "contour.hpp"
#include "tree.hpp"

namespace HMCont2D{ namespace Assembler{


//returns all possible contours (including 1 edge contours) from collection
//direction is arbitrary
//Resulting contours will not intersect each other if input edges don't 
std::vector<HMCont2D::Contour> AllContours(const HMCont2D::ECollection& input);

// =================== Extended Tree
//returns a tree which contains all edges from input
HMCont2D::ExtendedTree ETree(const HMCont2D::ECollection& input);


// =================== Single Contour
//Assemble single contour from shattered edges starting from given points of collection edges
HMCont2D::Contour Contour1(const ECollection& col, const Point* pnt_start, const Point* pnt_end);
HMCont2D::Contour Contour1(const ECollection& col, const Point* pnt_start);

//Assemble from another contour
HMCont2D::Contour Contour1(const Contour& col, const Point* pnt_start, const Point* pnt_end);
HMCont2D::Contour Contour1(const Contour& col, Point pnt_start, Point pnt_end);

//assemles for pnt_start in the direction (+-1) til the length of contour
//resulting contour will be longer or equal to givenn len
HMCont2D::Contour Contour1(const Contour& col, const Point* pnt_start, int direction, double len);

}}
#endif

