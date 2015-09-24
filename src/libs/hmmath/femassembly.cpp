#include "femassembly.hpp"

using namespace HMFem::Impl;

namespace{

// =========================== 3 node cells
shared_ptr<LocMat> LaplasLocalMatrix3(const vector<const GridPoint*>& points){
	assert(points.size() == 3);
	shared_ptr<LocMat3Sym> ret(new LocMat3Sym());
	const double &x1 = points[0]->x, &x2 = points[1]->x, &x3 = points[2]->x;
	const double &y1 = points[0]->y, &y2 = points[1]->y, &y3 = points[2]->y;
	double J11 = x2 - x1, J12 = x3 - x1;
	double J21 = y2 - y1, J22 = y3 - y1;
	double modJ = J22*J11 - J21*J12;
	double J2122 = J21 - J22, J1112 = J11 - J12;

	(*ret)[0] = ( sqr(J2122) + sqr(J1112));
	(*ret)[1] = ( J22*J2122  + J12*J1112 );
	(*ret)[2] = (-J11*J1112  - J21*J2122 );
	(*ret)[3] = ( sqr(J22)   + sqr(J12)  );
	(*ret)[4] = (-J21*J22    - J12*J11   );
	(*ret)[5] = ( sqr(J11)   + sqr(J21)  );

	return ret;
}

// ============================ 4 node cells
//integration in [-1, 1]x[-1, 1] square
struct _Integration{
	virtual double operator()(
		std::function<double(double, double)> fun) = 0;
};
//4nodes square integration
struct Integration4: public _Integration{

	const double x[4]={ 0.577350269189626e0,
			-0.577350269189626e0,
			-0.577350269189626e0,
			 0.577350269189626e0
	};
	const double y[4]={ 0.577350269189626e0,
			 0.577350269189626e0,
			-0.577350269189626e0,
			-0.577350269189626e0
	};
	const double w[4]={ 1,1,1,1};

	double operator()(std::function<double(double, double)> fun){
		_THROW_NOT_IMP_;
	}
};

//9nodes square integration
struct Integration9: public _Integration{
	const double x[9] = { 0.774596669241483e0,
			-0.774596669241483e0,
			 0.774596669241483e0,
			-0.774596669241483e0,
			 0.774596669241483e0,
			-0.774596669241483e0,
			 0.0e0,
			 0.0e0,
			 0.0e0
	};
	const double y[9] = { 
			 0.774596669241483e0,
			 0.774596669241483e0,
			-0.774596669241483e0,
			-0.774596669241483e0,
			 0.0e0,
			 0.0e0,
			 0.774596669241483e0,
			-0.774596669241483e0,
			 0.0e0
	};
	const double w[9] = {
			0.308641975308642e0,
			0.308641975308642e0,
			0.308641975308642e0,
			0.308641975308642e0,
			0.493827160493827e0,
			0.493827160493827e0,
			0.493827160493827e0,
			0.493827160493827e0,
			0.790123456790123e0
	};
	double operator()(std::function<double(double, double)> fun){
		return w[0] * fun(x[0], y[0]) +
		       w[1] * fun(x[1], y[1]) +
		       w[2] * fun(x[2], y[2]) +
		       w[3] * fun(x[3], y[3]) +
		       w[4] * fun(x[4], y[4]) +
		       w[5] * fun(x[5], y[5]) +
		       w[6] * fun(x[6], y[6]) +
		       w[7] * fun(x[7], y[7]) +
		       w[8] * fun(x[8], y[8]);
	}
};

//using 9-node integration
std::unique_ptr<_Integration> NumIntegration(new Integration9());

//stiff matrix 
double stiff4(const vector<const GridPoint*>& pts, int i, int j, double k, double e){
	//basic function derivatives
	static const std::function<double(double, double)> ddk[4] ={
		[](double x, double y){ return -(1-y)/4; },
		[](double x, double y){ return  (1-y)/4; },
		[](double x, double y){ return  (1+y)/4; },
		[](double x, double y){ return -(1+y)/4;}
	};
	static const std::function<double(double, double)> dde[4] ={
		[](double x, double y){ return -(1-x)/4; },
		[](double x, double y){ return -(1+x)/4; },
		[](double x, double y){ return  (1+x)/4; },
		[](double x, double y){ return  (1-x)/4; }
	};

	const double &x0 = pts[0]->x, &x1 = pts[1]->x, &x2 = pts[2]->x, &x3 = pts[3]->x;
	const double &y0 = pts[0]->y, &y1 = pts[1]->y, &y2 = pts[2]->y, &y3 = pts[3]->y;
	double j11 = ((1-e)*(x1-x0)+(1+e)*(x2-x3))/4.0;
	double j12 = ((1-k)*(x3-x0)+(1+k)*(x2-x1))/4.0;
	double j21 = ((1-e)*(y1-x0)+(1+e)*(y2-y3))/4.0;
	double j22 = ((1-k)*(y3-y0)+(1+k)*(y2-y1))/4.0;
	double modj = j11*j22-j12*j21;

	double didk = ddk[i](k, e); 
	double djdk = ddk[j](k, e); 
	double dide = dde[i](k, e); 
	double djde = dde[j](k, e); 

	double d1dx= j22*didk - j21*dide;
	double d1dy=-j12*didk + j11*dide;
	double d2dx= j22*djdk - j21*djde;
	double d2dy=-j12*djdk + j11*djde;

	return (d1dx*d2dx+d1dy*d2dy)/modj;
}

shared_ptr<LocMat> LaplasLocalMatrix4(const vector<const GridPoint*>& points){
	assert(points.size() == 4);
	using namespace std::placeholders;
	shared_ptr<LocMat4Sym> ret(new LocMat4Sym());
	auto& I = *NumIntegration;
	(*ret)[0] = I(std::bind(stiff4, points, 0, 0, _1, _2));
	(*ret)[1] = I(std::bind(stiff4, points, 0, 1, _1, _2));
	(*ret)[2] = I(std::bind(stiff4, points, 0, 2, _1, _2));
	(*ret)[3] = I(std::bind(stiff4, points, 0, 3, _1, _2));
	(*ret)[4] = I(std::bind(stiff4, points, 1, 1, _1, _2));
	(*ret)[5] = I(std::bind(stiff4, points, 1, 2, _1, _2));
	(*ret)[6] = I(std::bind(stiff4, points, 1, 3, _1, _2));
	(*ret)[7] = I(std::bind(stiff4, points, 2, 2, _1, _2));
	(*ret)[8] = I(std::bind(stiff4, points, 2, 3, _1, _2));
	(*ret)[9] = I(std::bind(stiff4, points, 3, 3, _1, _2));
	return ret;
}

shared_ptr<LocMat> LaplasLocalMatrix(const vector<const GridPoint*>& points){
	if (points.size() == 3) return LaplasLocalMatrix3(points);
	else return LaplasLocalMatrix4(points);
}


}

// ============================== Global assembling
Mat Assemble::PureLaplas(const Grid43& grid){
	Mat ret;
	for (int i=0; i<grid.n_cells(); ++i){
		vector<const GridPoint*> pp = grid.get_cell(i)->get_points();
		shared_ptr<LocMat> A = LaplasLocalMatrix(pp);
		vector<int> pind(pp.size());
		std::transform(pp.begin(), pp.end(), pind.begin(), 
				[](const GridPoint* p){ return p->get_ind(); });
		A->ToMat(pind, ret);
	}
	return ret;
}


