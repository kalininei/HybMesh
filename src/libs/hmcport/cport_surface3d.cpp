#include "hmproject.h"
#include "cport_surface3d.h"
#include "hmgrid3d.hpp"

using HMGrid3D::SSurface;

namespace{
HMGrid3D::SGrid& hmsurf(CPortSurface3D* g){
	return *static_cast<HMGrid3D::SSurface*>(g);
}
const HMGrid3D::SGrid& hmsurf(const CPortGrid3D* g){
	return *static_cast<const HMGrid3D::SSurface*>(g);
}
}

void* extract_grid3_surface(const void* grid3){
	SSruface* ret = 0;
	try{
		const HMGrid3D::SGrid* g = static_cast<const HMGrid3D::SGrid*>(grid3);
		auto sers = g->serialized_surface();
		ret = new SSurface(std::move(sers));
		return ret;
	} catch (const std::exception &e) {
		if (ret != 0) delete ret;
		std::cout<<e.what()<<std::endl;
		return NULL;
	}
}
