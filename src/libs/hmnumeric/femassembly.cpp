#include "femassembly.hpp"

using namespace HMFem;

namespace{

// =========================== 3 node cells
shared_ptr<HMMath::LocMat> LaplasLocalMatrix3(const HM2D::VertexData& points){
	assert(points.size() == 3);
	shared_ptr<HMMath::LocMat3Sym> ret(new HMMath::LocMat3Sym());
	const double &x1 = points[0]->x, &x2 = points[1]->x, &x3 = points[2]->x;
	const double &y1 = points[0]->y, &y2 = points[1]->y, &y3 = points[2]->y;

	// x(k, e) = k*(x2-x1) + e*(x3-x1) + x1
	// y(k, e) = k*(y2-y1) + e*(y3-y1) + y1
	//
	// J    = |dxdk dydk|
	//        |dxde dyde|
	//
	// InvJ =  1  | J22  -J12|
	//        |J| |-J21   J11|
	//
	// |dfdx| = InvJ * |dfdk|
	// |dfdy|          |dfde|

	double j11 = x2 - x1, j21 = x3 - x1;
	double j12 = y2 - y1, j22 = y3 - y1;
	double modjx2 = 2*(j22*j11 - j21*j12);
	double j1222 = j12 - j22, j1121 = j11 - j21;

	(*ret)[0] = ( sqr(j1222) + sqr(j1121))/modjx2;
	(*ret)[1] = ( j22*j1222  + j21*j1121 )/modjx2;
	(*ret)[2] = (-j11*j1121  - j12*j1222 )/modjx2;
	(*ret)[3] = ( sqr(j22)   + sqr(j21)  )/modjx2;
	(*ret)[4] = (-j12*j22    - j21*j11   )/modjx2;
	(*ret)[5] = ( sqr(j11)   + sqr(j12)  )/modjx2;

	return ret;
}

shared_ptr<HMMath::LocMat> DDxLocalMatrix3(const HM2D::VertexData& points){
	assert(points.size() == 3);
	shared_ptr<HMMath::LocMat3> ret(new HMMath::LocMat3());
	const double &y1 = points[0]->y, &y2 = points[1]->y, &y3 = points[2]->y;
	double j12 = y2 - y1, j22 = y3 - y1;

	(*ret)[0] = (*ret)[3] = (*ret)[6] = (-j22+j12)/6.0;
	(*ret)[1] = (*ret)[4] = (*ret)[7] = j22/6.0;
	(*ret)[2] = (*ret)[5] = (*ret)[8] = -j12/6.0;

	return ret;
}

shared_ptr<HMMath::LocMat> DDyLocalMatrix3(const HM2D::VertexData& points){
	assert(points.size() == 3);
	shared_ptr<HMMath::LocMat3> ret(new HMMath::LocMat3());
	const double &x1 = points[0]->x, &x2 = points[1]->x, &x3 = points[2]->x;

	double j11 = x2 - x1, j21 = x3 - x1;

	(*ret)[0] = (*ret)[3] = (*ret)[6] = (j21-j11)/6.0;
	(*ret)[1] = (*ret)[4] = (*ret)[7] = -j21/6.0;
	(*ret)[2] = (*ret)[5] = (*ret)[8] = j11/6.0;

	return ret;
}

shared_ptr<HMMath::LocMat> FullMassLocalMatrix3(const HM2D::VertexData& points){
	assert(points.size() == 3);
	shared_ptr<HMMath::LocMat3Sym> ret(new HMMath::LocMat3Sym());
	const double &x1 = points[0]->x, &x2 = points[1]->x, &x3 = points[2]->x;
	const double &y1 = points[0]->y, &y2 = points[1]->y, &y3 = points[2]->y;
	double j11 = x2 - x1, j21 = x3 - x1;
	double j12 = y2 - y1, j22 = y3 - y1;
	double v1 = (j22*j11 - j21*j12)/12.0;
	double v2 = v1/2.0;

	(*ret)[0] = v1;
	(*ret)[1] = v2;
	(*ret)[2] = v2;
	(*ret)[3] = v1;
	(*ret)[4] = v2;
	(*ret)[5] = v1;

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
double stiff4(const HM2D::VertexData& pts, int i, int j, double k, double e){
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
	double j21 = ((1-k)*(x3-x0)+(1+k)*(x2-x1))/4.0;
	double j12 = ((1-e)*(y1-y0)+(1+e)*(y2-y3))/4.0;
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

shared_ptr<HMMath::LocMat> LaplasLocalMatrix4(const HM2D::VertexData& points){
	using namespace std::placeholders;
	assert(points.size() == 4);
	shared_ptr<HMMath::LocMat4Sym> ret(new HMMath::LocMat4Sym());
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

shared_ptr<HMMath::LocMat> DDxLocalMatrix4(const HM2D::VertexData& points){
	_THROW_NOT_IMP_;
}
shared_ptr<HMMath::LocMat> DDyLocalMatrix4(const HM2D::VertexData& points){
	_THROW_NOT_IMP_;
}
shared_ptr<HMMath::LocMat> FullMassLocalMatrix4(const HM2D::VertexData& points){
	_THROW_NOT_IMP_;
}

shared_ptr<HMMath::LocMat> LaplasLocalMatrix(const HM2D::VertexData& points){
	if (points.size() == 3) return LaplasLocalMatrix3(points);
	else return LaplasLocalMatrix4(points);
}
shared_ptr<HMMath::LocMat> FullMassLocalMatrix(const HM2D::VertexData& points){
	if (points.size() == 3) return FullMassLocalMatrix3(points);
	else return FullMassLocalMatrix4(points);
}
shared_ptr<HMMath::LocMat> DDxLocalMatrix(const HM2D::VertexData& points){
	if (points.size() == 3) return DDxLocalMatrix3(points);
	else return DDxLocalMatrix4(points);
}
shared_ptr<HMMath::LocMat> DDyLocalMatrix(const HM2D::VertexData& points){
	if (points.size() == 3) return DDyLocalMatrix3(points);
	else return DDyLocalMatrix4(points);
}


}

// ============================== Global assembling
namespace{

shared_ptr<HMMath::Mat> GlobAssembly(const HM2D::GridData& grid, 
		decltype(LaplasLocalMatrix)& fun){
	aa::enumerate_ids_pvec(grid.vvert);
	shared_ptr<HMMath::Mat> ret(new HMMath::Mat());
	for (int i=0; i<grid.vcells.size(); ++i){
		HM2D::VertexData pp = HM2D::Contour::OrderedPoints1(grid.vcells[i]->edges);
		shared_ptr<HMMath::LocMat> A = fun(pp);
		vector<int> pind(pp.size());
		std::transform(pp.begin(), pp.end(), pind.begin(), 
				[](const shared_ptr<HM2D::Vertex>& p){ return p->id; });
		A->ToMat(pind, *ret);
	}
	return ret;
}

}
shared_ptr<HMMath::Mat> Assemble::PureLaplas(const HM2D::GridData& grid){
	return GlobAssembly(grid, LaplasLocalMatrix);
}

shared_ptr<HMMath::Mat> Assemble::FullMass(const HM2D::GridData& grid){
	return GlobAssembly(grid, FullMassLocalMatrix);
}

shared_ptr<HMMath::Mat> Assemble::DDx(const HM2D::GridData& grid){
	return GlobAssembly(grid, DDxLocalMatrix);
}

shared_ptr<HMMath::Mat> Assemble::DDy(const HM2D::GridData& grid){
	return GlobAssembly(grid, DDyLocalMatrix);
}

vector<double> Assemble::LumpMass(const HM2D::GridData& grid){
	shared_ptr<HMMath::Mat> m = FullMass(grid);
	vector<double> tmp(m->rows(), 1.0);
	vector<double> ret(m->rows(), 0.0);
	m->MultVec(tmp, ret);
	return ret;
}

