#ifndef HMCONT2D_CONT_ASSEMBLER_HPP
#define HMCONT2D_CONT_ASSEMBLER_HPP

#include "contour.hpp"

namespace HMCont2D{ namespace Assembler{


//returns all possible contours (including 1 edge contours) from collection
std::vector<HMCont2D::Contour> AllContours(const HMCont2D::ECollection& input);


}}
#endif

