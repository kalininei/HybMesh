#include "tscaler.hpp"

using namespace Autoscale;

D2::D2(vector<Point>& p, double a){
	sc = ScaleBase::doscale(p.begin(), p.end(), a);
}
//constructors with automatic unscaling on destruction
D2::D2(HM2D::EdgeData* vd, double a){
	_einp.push_back(vd);
	auto av = HM2D::AllVertices(*vd);
	sc = HM2D::Scale01(av, a);
}
D2::D2(HM2D::GridData* vd, double a){
	_ginp.push_back(vd);
	sc = HM2D::Scale01(vd->vvert, a);
}
D2::D2(const vector<HM2D::EdgeData*>& vd, double a){
	if (vd.size() == 0) return;
	HM2D::VertexData av;
	for (auto it: vd){
		auto av1 = HM2D::AllVertices(*it);
		av.insert(av.end(), av1.begin(), av1.end());
		_einp.push_back(it);
	}
	sc = HM2D::Scale01(av, a);
}

D2::~D2(){
	for (auto it: _ginp){ HM2D::Unscale(it->vvert, sc); }
	for (auto it: _einp){ HM2D::Unscale(*it, sc); }
}

void D2::add_data(HM2D::EdgeData* vd){
	_einp.push_back(vd);
	HM2D::Scale(*vd, sc);
}
void D2::add_data(HM2D::GridData* vd){
	_ginp.push_back(vd);
	HM2D::Scale(vd->vvert, sc);
}
void D2::add_data(const vector<HM2D::EdgeData*>& vd){
	_einp.insert(_einp.end(), vd.begin(), vd.end());
	for (auto it: vd){ HM2D::Scale(*it, sc); }
}

void D2::scale(vector<Point>& pts){
	sc.scale(pts.begin(), pts.end());
}
void D2::scale(vector<double>& lens){
	for (auto& v: lens){ v/=sc.L; }
}
void D2::scale(Point& p){
	sc.scale(p);
}
void D2::scale(double& a){
	a/=sc.L;
}

void D2::unscale(HM2D::EdgeData* vd){
	HM2D::Unscale(*vd, sc);
}
void D2::unscale(HM2D::GridData* vd){
	HM2D::Unscale(vd->vvert, sc);
}
void D2::unscale(HM3D::GridData* vd){
	for (auto v: vd->vvert){
		v->z *= sc.L;
		v->x *= sc.L; v->x += sc.p0.x;
		v->y *= sc.L; v->y += sc.p0.y;
	}
}

D3::~D3(){
	for (auto it: _ginp) HM3D::Unscale(it->vvert, sc);
	for (auto it: _finp) HM3D::Unscale(*it, sc);
}
D3::D3(HM3D::FaceData* vd, double a){
	sc = HM3D::Scale01(*vd, a);
	_finp.push_back(vd);
}
D3::D3(HM3D::GridData* vd, double a){
	sc = HM3D::Scale01(vd->vvert, a);
	_ginp.push_back(vd);
}
D3::D3(const vector<HM3D::GridData*>& vd, double a){
	HM3D::VertexData av;
	for (auto it: vd){
		av.insert(av.end(), it->vvert.begin(), it->vvert.end());
	}
	sc = HM3D::Scale01(av, a);
	_ginp = vd;
}

D3::D3(const vector<HM3D::FaceData*>& vd, double a){
	HM3D::VertexData av;
	for (auto it: vd){
		HM3D::VertexData av1 = HM3D::AllVertices(*it);
		av.insert(av.end(), av1.begin(), av1.end());
	}
	sc = HM3D::Scale01(av, a);
	_finp = vd;
}

void D3::unscale(HM3D::GridData* g){
	HM3D::Unscale(g->vvert, sc);
}
