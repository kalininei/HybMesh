#include "finder3d.hpp"


std::tuple<int, double>
HM3D::Finder::ClosestPoint(const HM3D::VertexData& vec, const Point3& v){
	std::tuple<int, double> ret;
	int& ind = std::get<0>(ret);
	double& meas = std::get<1>(ret);
	if (vec.size() == 0){
		ind = -1; meas = -1;
		return ret;
	}

	ind = 0; meas = Point3::meas(*vec[0], v);
	for (int i=1; i<vec.size(); ++i){
		double m = Point3::meas(v, *vec[i]);
		if (m<meas){
			meas = m; ind = i;
		}
	}

	return ret;
}

