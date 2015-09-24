#include "hmcompute.hpp"
using namespace HMMath::Compute;

double HMMath::Compute::SecantMethod(double a, double b, const std::function<double(double)>& fun,
		int MaxIt, double eps){
	if (a<b) std::swap(a, b);
	else if (a==b) throw SecantMethodError("Secant method error: a == b");

	double fa = fun(a);
	if (fabs(fa) < eps) return a;
	double fb;

	int it = 0;
	while(fabs(b - a) > eps){
		fb = fun(b);
		a = b - (b-a)*fb/(fb-fa);
		fa = fun(a);
		b = a - (a-b)*fa/(fa-fb);

		if (++it > MaxIt) throw
			SecantMethodError("Secant method error: max iter was reached");
	}

	return (a+b)/2.0;
}


double HMMath::Compute::GDescentMin(double x0, const std::function<double(double)>& fun,
		double hgrad, double lambda,
		int MaxIt, double eps){
	int it = 0;
	double grad;
	double nrm, diff;

	while (it<MaxIt){
		++it;
		//find gradient
		grad = (fun(x0 + hgrad) - fun(x0 - hgrad))/(2*hgrad);
		//make a step
		diff = lambda*grad;
		x0 -= diff;
		//exit conditions
		if (fabs(diff) < eps) break;
	}

	if (it>=MaxIt) throw GDescentError("Gradient descent failed to converge");

	return x0;
}


namespace{

const double GoldRat=(sqrt(5.0)+1.0)/2.0;
const double GoldRatI=1.0/GoldRat;

}

double HMMath::Compute::GoldenRatioMin(double x1, double x2, const std::function<double(double)>& fun,
		int MaxIt, double eps){
	int it = 0;
	double x, x1_new, x2_new, f_x1, f_x2;

	while (it < MaxIt){
		x=(x1+x2)/2.0;
		double dif=fabs(x2-x1);
		if (dif<eps) break;  
		double v=dif*GoldRatI;
		x1_new=x2-v; x2_new=x1+v;
		f_x1=fun(x1_new);
		f_x2=fun(x2_new);
		(f_x1>=f_x2) ? x1=x1_new : x2=x2_new;
		++it;
	}

	//if (it>=MaxIt) throw GoldenRatioError("Golden ratio minimization failed to converge");
	return x;
}
