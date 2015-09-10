#include "simple_bgrid.hpp"

using namespace HMBlay::Impl;

SimpleBGrid::SimpleBGrid(ExtPath& pth, int _Nsmooth):BGrid(), 
		Nsmooth(_Nsmooth), _is_closed(pth.is_closed()), _hor_size(pth.size()){
	//Only NO and REGULAR corners in pth
	assert(pth.ext_data.size() == pth.ordered_points().size());
	assert(std::all_of(pth.ext_data.begin(), pth.ext_data.end(),
		[](PathPntData& d){
			return (d.tp == CornerTp::NO ||
				d.tp == CornerTp::REGULAR);
		}
	));
	// creates rectangular stencil
	FillStencil(pth);
	//assemble points
	points = stencil.GetAll();
	// assemble cells
	AssembleFromStencil();
	//post
	set_indicies();
}

void SimpleBGrid::FillStencil(ExtPath& pth){
	//build along normals
	stencil.Fill(pth);
	//guarantee no intersections
	stencil.NoIntersections();
	//smoothing
	for (int i=0; i<Nsmooth; ++i) stencil.SmoothStep(Wsmooth);
}

void SimpleBGrid::AssembleFromStencil(){
	for (int edind = 0; edind < hor_size(); ++edind){
		int stenind1 = edind;
		int stenind2 = (is_closed() && edind == hor_size() - 1) ? 0 : edind+1;
		//number of cells for this edge
		int n_height = std::min(stencil[stenind1].size(), stencil[stenind2].size()) - 1;
		for (int ih = 0; ih < n_height; ++ih){
			Cell* c = aa::add_shared(cells, Cell());
			add_point_to_cell(c, stencil[stenind1][ih].get());
			add_point_to_cell(c, stencil[stenind2][ih].get());
			add_point_to_cell(c, stencil[stenind2][ih+1].get());
			add_point_to_cell(c, stencil[stenind1][ih+1].get());
			priority[c] = HIGHEST_PRIORITY - ih;
		}
	}
}


void SimpleBGrid::TStencil::Fill(ExtPath& pth){
	is_closed = pth.is_closed();
	auto basepnts = pth.ordered_points();
	for (int i=0; i<pth.size() + (is_closed?0:1); ++i){
		Point* b = basepnts[i];
		Vect n = pth.ext_data[i].normal;
		auto& part = pth.ext_data[i].opt->partition;
		this->push_back(ShpVector<GridPoint>());
		for (int j = 0; j<part.size(); ++j){
			aa::add_shared(this->back(), GridPoint(*b + n * part[j]));
		}
	}
	fill_angles();
	fill_basedist();
}

ShpVector<GridPoint> SimpleBGrid::TStencil::GetAll(){
	ShpVector<GridPoint> ret;
	for (auto i=0; i<size(); ++i){
		auto& e = operator[](i);
		std::copy(e.begin(), e.end(), std::back_inserter(ret));
	}
	return ret;
}

void SimpleBGrid::TStencil::NoIntersections(){
	auto has_intersections = [&]()->bool{
		const double crosseps = 0.01;
		double ksieta[2];
		for (auto it=begin(); it!=end(); ++it){
			iterator itnext = std::next(it);
			if (it == end() - 1){
				if (is_closed) itnext = begin();
				else break;
			}
			Point *p0 = (*it)[0].get(), *p1 = (*it).back().get();
			Point *p2 = (*itnext)[0].get(), *p3 = (*itnext).back().get();
			SectCross(*p0, *p1, *p2, *p3, ksieta);
			if (ksieta[0]<1+crosseps && ksieta[0]>-crosseps) return true;
			if (ksieta[1]<1+crosseps && ksieta[1]>-crosseps) return true;
		}
		return false;
	};

	int it = 0;
	int itmax = 10;
	while(has_intersections() && it<itmax) {SmoothStep(0.0); it++;}

	if (it==itmax){
		std::cout<<"Failed to remove normal intersections in a regular region"<<std::endl;
		assert(false);
	}
}

void SimpleBGrid::TStencil::SmoothStep(double w){
	for (int i=0; i<size(); ++i){
		if (!is_closed && (i==0 || i==size()-1)) continue;
		int inext = i+1, iprev = i-1;
		if (is_closed && inext==size()) inext = 0;
		if (is_closed && iprev<0) iprev = size()-1;
		
		double angle_cur = get_angle(i);
		double angle_prev = get_angle(iprev);
		double angle_next = get_angle(inext);
		//continue if angles are equal
		if (fabs(angle_cur-angle_prev)<1e-3 &&
			fabs(angle_cur-angle_next)<1e-3) continue;

		double rot1 = ToAngle(angle_prev-angle_cur);
		double rot2 = ToAngle(angle_next-angle_cur);
		if (rot1 > M_PI) rot1-=2*M_PI;
		if (rot2 > M_PI) rot2-=2*M_PI;
		double k1=basedist(i, iprev), k2=basedist(i, inext);
		k1 = k1/(k1+k2); k2=1.0-k1;
		k1=0.5; k2=0.5;
		double rot = (1.0-w)*(k2*rot2+k1*rot1);

		set_angle(i, angle_cur + rot);
	}
}

void SimpleBGrid::TStencil::set_angle(int i, double a){
	a = ToAngle(a);
	double rotangle = a - get_angle(i);
	auto& e = operator[](i);
	for (int k=1; k<e.size(); ++k){
		auto v = vecRotate(*e[k]-*e[0], rotangle);
		e[k]->x = e[0]->x+v.x;
		e[k]->y = e[0]->y+v.y;
	}
	_angles[i] = a;
}

void SimpleBGrid::TStencil::fill_angles(){
	_angles.resize(size());
	for (int i=0; i<size(); ++i){
		auto& e = operator[](i);
		_angles[i] = ToAngle(atan2(e[1]->y-e[0]->y, e[1]->x-e[0]->x));
	}
}

void SimpleBGrid::TStencil::fill_basedist(){
	_basedist.resize(size());
	for (int i=0; i<size(); ++i){
		int inext = i+1;
		if (inext == size()){
			if (!is_closed) break;
			else inext = 0;
		}
		_basedist[i] = Point::dist(*this->operator[](i)[0], *this->operator[](inext)[0]);
	}
}

double SimpleBGrid::TStencil::basedist(int i, int inext) const{
	assert(abs(i-inext)==1 || abs(i-inext)==size()-1);
	int ia = std::min(i, inext);
	return _basedist[ia];
}




