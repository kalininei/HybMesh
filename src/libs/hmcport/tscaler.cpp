#include "tscaler.hpp"

using namespace Autoscale;

void D3::_add_to_cache(const HM3D::VertexData& vd){
	_cached.insert(_cached.end(), vd.begin(), vd.end());
}

D3::D3(const HM3D::VertexData& vd, double a): side(a){
	_cached = vd;
	_do_scale();
}
D3::D3(const HM3D::FaceData& fd, double a): side(a){
	auto av = AllVertices(fd);
	_cached=av;
	_do_scale();
}
D3::D3(std::initializer_list<const HM3D::VertexData*> vlst, double a): side(a){
	for (auto v: vlst) _add_to_cache(*v);
	_do_scale();
}

D3::~D3(){
	_undo_scale();
}
void D3::add_data(HM3D::VertexData& vd){
	scale(vd);
	_add_to_cache(vd);
}
void D3::add_data(HM3D::FaceData& vd){
	auto av = AllVertices(vd);
	scale(av);
	_add_to_cache(av);
}

void D3::_do_scale(){
	sc = ScaleBase3::doscale(_cached, side);
}

void D3::_undo_scale(){
	sc.unscale(_cached.begin(), _cached.end());
}

void D3::scale(vector<double>& vd){
	for (auto& a: vd) a/=sc.L;
}
void D3::scale(HM3D::VertexData& vd){
	sc.scale(vd.begin(), vd.end());
}
void D3::unscale(HM3D::VertexData& vd){
	sc.unscale(vd.begin(), vd.end());
}
