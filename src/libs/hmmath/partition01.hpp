#ifndef HYBMESH_MATH_PARTITION01_HPP
#define HYBMESH_MATH_PARTITION01_HPP
#include "hmproject.h"

namespace HMMath{

// from [length -> recommended step size] map
// to discretized [0, 1] section.
vector<double> Partition01(const std::map<double, double>& wts);

//rounds floating number vector values keeping its sum.
vector<int> RoundVector(const vector<double>& vect, const vector<int>& minsize);

}

#endif
