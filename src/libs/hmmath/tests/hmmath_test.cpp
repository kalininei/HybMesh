#include "piecewise.hpp"
#include "hmtesting.hpp"
using HMTesting::add_check;

void test01(){
	std::cout<<"01. Piecewise testing"<<std::endl;
	HMMath::StepFunction<double> f1(0);
	f1.add_section(-2, -1, 7.2);
	add_check(ISEQ(f1.Integral(-3, 0), 7.2) &&
		ISEQ(f1.Integral(-2, -1), 7.2) &&
		ISEQ(f1.Integral(-2, -1.9), 0.72) &&
		ISEQ(f1.Integral(-1.6, -1.5), 0.72) &&
		ISEQ(f1.Integral(-1.2, -1), 1.44) &&
		ISEQ(f1.Integral(3, 5), 0), "step function 1");
	f1.add_section(-1, 0, 3);
	f1.add_section(-1, -0.55, 1);
	f1.add_section(-1.8, 5, 12);
	f1.add_section(-10, 0, 22.1);
	add_check(ISEQ(f1.Integral(-5, 5), 170.5), "step function 2");
	f1.add_section(-f1.inf(), f1.inf(), 1);
	f1.add_section(-1, 1, 7);
	f1.add_section(-1+1e-12, -0.5, 8);
	add_check(ISEQ(f1.Integral(-10, 10), 32.5), "step function 3");
	add_check(ISEQ(f1(-10), 1) &&
		ISEQ(f1(-1), 8) &&
		ISEQ(f1(0), 7) &&
		ISEQ(f1(0.5), 7) &&
		ISEQ(f1(100), 1), "step function 4");

	HMMath::LinearPiecewise f2;
	f2.add_point(1, 1);
	f2.add_point(0.5, 0.5);
	f2.add_point(0, 0);
	f2.add_point(2, 0);
	add_check(ISEQ(f2(-1), -1) && ISEQ(f2(0), 0) && ISEQ(f2(0.3), 0.3) && ISEQ(f2(1), 1)
			&& ISEQ(f2(10), -8), "linear piecewise 1");
	add_check(ISEQ(f2.Integral(-10, 10), -81) &&
		ISEQ(f2.Integral(-1, 1), 0), "linear piecewise 2");
}

int main(){
	test01();

	HMTesting::check_final_report();
	std::cout<<"DONE"<<std::endl;
}
