#include "hmproject.h"
#include "cport_surface3d.h"
#include "hmgrid3d.hpp"

using HMGrid3D::SSurface;

namespace{
HMGrid3D::SSurface& hmsurf(CPortSurface3D* g){
	return *static_cast<HMGrid3D::SSurface*>(g);
}
const HMGrid3D::SSurface& hmsurf(const CPortSurface3D* g){
	return *static_cast<const HMGrid3D::SSurface*>(g);
}
}

void free_srf3(CPortSurface3D* s){
	delete static_cast<HMGrid3D::SSurface*>(s);
}

void surf3_dims(const CPortSurface3D* s, int* dims){
	dims[0] = hmsurf(s).n_vert;
	dims[1] = hmsurf(s).n_edges;
	dims[2] = hmsurf(s).n_faces;
}

int surf3_btypes(const CPortSurface3D* s3, int* ret){
	try{
		auto s = hmsurf(s3);
		std::copy(s.btypes.begin(), s.btypes.begin() + s.n_faces, ret);
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return 0;
	}
	return 1;
}

int surf3_volume(const CPortSurface3D* s, double* ret){
	try{
		auto& s3 = hmsurf(s);
		auto tree = HMGrid3D::SurfaceTree::Assemble(s3);
		auto rv = HMGrid3D::SurfTreeTReverter(tree);
		*ret = fabs(HMGrid3D::Surface::Volume(s3));
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return 0;
	}
	return 1;
}

void* extract_grid3_surface(const void* grid3){
	SSurface* ret = 0;
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

int extract_subsurfaces(const CPortSurface3D* src, int* num, CPortSurface3D*** ret){
	try{
		auto srf = hmsurf(src);
		vector<HMGrid3D::Surface> subs = HMGrid3D::Surface::AllSubSurfaces(srf);
		*num = subs.size();
		*ret = new void*[*num];
		for (int i=0; i<*num; ++i){
			(*ret)[i] = new SSurface(std::move(subs[i]));
		}
	} catch (const std::exception &e) {
		std::cout<<e.what()<<std::endl;
		return 0;
	}
	return 1;
}

void* s3writer_create(const char* sname, CPortSurface3D* surf, void* awriter,
		void* subnode, const char* fmt){
	try{
		HMXML::ReaderA* wr = static_cast<HMXML::ReaderA*>(awriter);
		HMXML::Reader* sn = static_cast<HMXML::Reader*>(subnode);
		SSurface* s = static_cast<SSurface*>(surf);
		auto ret = new HMGrid3D::Export::SurfaceWriter(*s, wr, sn, sname, fmt);
		return ret;
	} catch (const std::exception &e){
		std::cout<<e.what()<<std::endl;
		return 0;
	}
}
int s3writer_add_defined_field(void* swriter, const char* field){
	try{
		_THROW_NOT_IMP_;
	} catch (const std::exception &e) {
		std::cout<<e.what()<<std::endl;
		return 0;
	}
	return 1;
}
void s3writer_free(void* swriter){
	try{
		delete static_cast<HMGrid3D::Export::SurfaceWriter*>(swriter);
	} catch (const std::exception &e){
		std::cout<<e.what()<<std::endl;
	}
}

void* s3reader_create(void* awriter, void* subnode, char* outname){
	try{
		auto wr = static_cast<HMXML::ReaderA*>(awriter);
		auto sn = static_cast<HMXML::Reader*>(subnode);
		//name
		std::string nm = sn->attribute(".", "name");
		if (nm.size()>1000) throw std::runtime_error("surface name is too long: " + nm);
		strcpy(outname, nm.c_str());
		//reader
		auto ret = HMGrid3D::Import::ReadHMC(wr, sn);
		return ret.release();
	} catch (const std::exception &e){
		std::cout<<e.what()<<std::endl;
		return 0;
	}
}
CPortSurface3D* s3reader_getresult(void* rd){
	try{
		auto reader = static_cast<HMGrid3D::Import::SurfaceReader*>(rd);
		return reader->result.release();
	} catch (const std::exception &e){
		std::cout<<e.what()<<std::endl;
		return 0;
	}
}
void s3reader_free(void* sreader){
	try{
		delete (HMGrid3D::Import::SurfaceReader*)sreader;
	} catch (const std::exception &e){
		std::cout<<e.what()<<std::endl;
	}
}
