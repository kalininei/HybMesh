#include "laplace_bem2d.hpp"
#include "densemat.hpp"
#include "treverter2d.hpp"
#include "hmtimer.hpp"

using namespace HMBem;

LaplaceCE2D::LaplaceCE2D(const HM2D::Contour::Tree& area){
	no_bmap = true;
	N = 0;
	for (auto& nd: area.nodes){
		//implemented only for single root - mutlitple children trees without cuts
		assert(area.roots().size() == 1 && nd->isbound() && nd->level<2);
		N+=nd->contour.size();
	}
	iprev.resize(N);
	dirvals.resize(N, 0);
	neuvals = len = len2 = loglen = nx = ny = xs = ys = xm = ym = dirvals;
	isdir.resize(N, false);
	int it=0;
	for (auto& nd: area.nodes){
		HM2D::Contour::R::ReallyDirect rev(nd->contour);
		bool counterclockwise = (HM2D::Contour::Area(nd->contour) > 0);
		bool correct_dir = nd->isouter() ? counterclockwise : !counterclockwise;
		for (auto& e: nd->contour){
			iprev[it] = it-1;
			len[it] = e->length();
			len2[it] = len[it]*len[it];
			loglen[it] = log(len[it]);
			xm[it] = (e->pfirst()->x + e->plast()->x)/2.0;
			ym[it] = (e->pfirst()->y + e->plast()->y)/2.0;
			if (correct_dir){
				xs[it] = e->pfirst()->x;
				ys[it] = e->pfirst()->y;
				nx[it] = (e->plast()->y - ys[it])/len[it];
				ny[it] = (xs[it] - e->plast()->x)/len[it];
			} else {
				xs[it] = e->plast()->x;
				ys[it] = e->plast()->y;
				nx[it] = (e->pfirst()->y - ys[it])/len[it];
				ny[it] = (xs[it] - e->pfirst()->x)/len[it];
			}
			++it;
		}
		iprev[it-nd->contour.size()] = it-1;
	}
};

void LaplaceCE2D::Solve(){
	HMMath::DenseMat mat(N);
	std::vector<double> rhs(N, 0);
	vector<double>::iterator matiter = mat.dt.begin();
	for (int i=0; i<N; ++i)
	for (int j=0; j<N; ++j){
		double dirac = (i==j) ? 0.5 : 0.0;
		auto F = Ffun(j, xm[i], ym[i]);
		if (isdir[j]){
			*matiter++ = -F.first;
			rhs[i] += dirvals[j]*(-F.second+dirac);
		} else{
			*matiter++ = F.second - dirac;
			rhs[i] += neuvals[j]*F.first;
		}
	}

	HMMath::DenseSolver slv(mat);
	std::vector<double> z(N, 0);
	slv.Solve(rhs, z);

	for (int i=0; i<N; ++i){
		if (isdir[i]) neuvals[i] = z[i];
		else dirvals[i] = z[i];
	}
}

int LaplaceCE2D::where_is(double x, double y) const{
	_THROW_NOT_IMP_;
}

double LaplaceCE2D::value_at(double x, double y) const{
	int w = where_is(x, y);
	if (w == INSIDE) return boundary_value_at(x, y);
	else if (w==BOUND) return internal_value_at(x, y);
	else throw std::runtime_error("point is outside given domain");
}

double LaplaceCE2D::internal_value_at(double x, double y) const{
	double ret = 0;
	for (int i=0; i<N; ++i){
		auto v = Ffun(i, x, y);
		ret += dirvals[i]*v.second-neuvals[i]*v.first;
	}
	return ret;
}

double LaplaceCE2D::boundary_value_at(double x, double y) const{
	if (no_bmap){
		for (int i=0; i<N; ++i) vfinder.add(xs[i], ys[i], i);
		no_bmap = false;
	}
	double ret = 0;
	for (int i=0; i<N; ++i){
		auto v = Ffun(i, x, y);
		ret += dirvals[i]*v.second-neuvals[i]*v.first;
	}

	auto fnd = vfinder.find(x, y);
	if (fnd != vfinder.end()){
		int i2 = fnd.data();
		int i1 = iprev[i2];
		double phi = Angle(Point(xm[i1], ym[i1]), Point(xs[i2], ys[i2]), Point(xm[i2], ym[i2]));
		return ret*2.0*M_PI/phi;
	} else {
		return ret/2.0;
	}
	return ret;
}

std::pair<double, double> LaplaceCE2D::Ffun(int k, double x, double y) const{
	double b = B(k, x, y);
	double e = E(k, x, y);
	double D = 4*len2[k]*e-b*b;
	double ba = b/2.0/len2[k];
	if (D < 1e-12){
		double balog1 = (ba == -1.0) ? 0 : (ba+1.0)*log(fabs(ba+1.0));
		double balog2 = (ba == 0.0) ? 0 : ba*log(fabs(ba));
		return std::make_pair(
			len[k]/2.0/M_PI*(loglen[k]+balog1-balog2-1.0),
			0.0);
	} else {
		double sq = sqrt(D);
		double at1 = atan((2.0*len2[k]+b)/sq) - atan(b/sq);
		return std::make_pair(
			len[k]/4.0/M_PI*(2.0*(loglen[k]-1.0) - ba*log(fabs(e/len2[k]))+(1.0+ba)*log(fabs(1.0+2.0*ba+e/len2[k]))+sq/len2[k]*at1),
			len[k]*(nx[k]*(xs[k]-x)+ny[k]*(ys[k]-y))/M_PI/sq*at1);
	}

}
double LaplaceCE2D::B(int k, double x, double y) const{
	return (-ny[k]*(xs[k]-x)+nx[k]*(ys[k]-y))*2.0*len[k];
}
double LaplaceCE2D::E(int k, double x, double y) const{
	return (xs[k]-x)*(xs[k]-x)+(ys[k]-y)*(ys[k]-y);
}
