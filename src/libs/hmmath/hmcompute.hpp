#ifndef HYBMESH_HMMATH_HMCOMPUTE_HPP
#define HYBMESH_HMMATH_HMCOMPUTE_HPP
#include "hmproject.h"

namespace HMMath{ namespace Compute{

// ======================= Secant Method =========================================
//throws if MaxIt was reached
class SecantMethodError: public std::runtime_error{
public:
	SecantMethodError(const char* s) noexcept: std::runtime_error(s){}
};
double SecantMethod(double x0, double x1, const std::function<double(double)>& fun,
		int MaxIt, double eps);


// ======================= Gradient descent minimization of 1d function
class GDescentError: public std::runtime_error{
public:
	GDescentError(const char* s) noexcept: std::runtime_error(s){}
};
//gradient descent with unknown gradient function
double GDescentMin(double x0, const std::function<double(double)>& fun,
		double hgrad, double lambda,
		int MaxIt, double eps);


// ======================== Golden ratio minimization
class GoldenRatioError: public std::runtime_error{
public:
	GoldenRatioError(const char* s) noexcept: std::runtime_error(s){}
};
double GoldenRatioMin(double a, double b, const std::function<double(double)>& fun,
		int MaxIt, double eps);



}}



#endif
