#ifndef HYBMESH_HMMATH_HMCOMPUTE_HPP
#define HYBMESH_HMMATH_HMCOMPUTE_HPP
#include "hmproject.h"

namespace HMMath{ namespace Compute{

// ======================= simple matrix operations
inline double determinant_2x2(double* A){
	return A[0]*A[3]-A[1]*A[2];
}

inline bool mat_inverse_2x2(double* A, double* B){
	double det=determinant_2x2(A);
	if (ISZERO(det)) return false;
	B[0]= A[3]/det;
	B[1]=-A[1]/det;
	B[2]=-A[2]/det;
	B[3]= A[0]/det;
	return true;
}

inline void mat_mult_vec_2x2(double* A, double* x, double* b){
	b[0] = A[0]*x[0] + A[1]*x[1];
	b[1] = A[2]*x[0] + A[3]*x[1];
}

inline double determinant_3x3(double* A){
	return A[0]*(A[4]*A[8]-A[5]*A[7])+
	       A[1]*(A[5]*A[6]-A[3]*A[8])+
	       A[2]*(A[3]*A[7]-A[4]*A[6]);
};

inline void mat_mult_vec_3x3(double* A, double* x, double* b){
	b[0] = A[0]*x[0] + A[1]*x[1] + A[2]*x[2];
	b[1] = A[3]*x[0] + A[4]*x[1] + A[5]*x[2];
	b[2] = A[6]*x[0] + A[7]*x[1] + A[8]*x[2];
}

inline bool mat_inverse_3x3(double* A, double* B){
	double det=determinant_3x3(A);
	if (ISZERO(det)) return false;
	B[0] = (A[4]*A[8]-A[5]*A[7])/det;
	B[1] = (A[2]*A[7]-A[1]*A[8])/det;
	B[2] = (A[1]*A[5]-A[2]*A[4])/det;
	B[3] = (A[5]*A[6]-A[3]*A[8])/det;
	B[4] = (A[0]*A[8]-A[2]*A[6])/det;
	B[5] = (A[2]*A[3]-A[0]*A[5])/det;
	B[6] = (A[3]*A[7]-A[4]*A[6])/det;
	B[7] = (A[1]*A[6]-A[0]*A[7])/det;
	B[8] = (A[0]*A[4]-A[1]*A[3])/det;
	return true;
}


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
