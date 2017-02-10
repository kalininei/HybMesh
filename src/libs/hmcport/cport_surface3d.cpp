#include "hmcport.h"
#include "cport_surface3d.h"
#include "c2cpp_helper.hpp"
#include "primitives3d.hpp"
#include "surface_tree.hpp"
#include "treverter3d.hpp"
#include "export3d_hm.hpp"

int s3_dims(void* obj, int* dims){
	try{
		auto surf = static_cast<HM3D::FaceData*>(obj);
		for (auto& f: *surf){
			aa::constant_ids_pvec(f->edges, 0);
			for (auto& e: f->edges){
				e->first()->id = 0;
				e->last()->id = 0;
			}
		}
		dims[0] = 0;
		dims[1] = 0;
		dims[2] = surf->size();
		for (auto& f: *surf)
		for (auto& e: f->edges) if (e->id == 0){
			++dims[1];
			e->id = 1;
			if (e->first()->id == 0){
				++dims[0];
				e->first()->id = 1;
			}
			if (e->last()->id == 0){
				++dims[0];
				e->last()->id = 1;
			}
		}
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}
int s3_tab_btypes(void* obj, int* ret){
	try{
		auto surf = static_cast<HM3D::FaceData*>(obj);
		for (auto& f: *surf) *ret++ = f->boundary_type;
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}
int s3_area(void* obj, double* ret){
	try{
		_THROW_NOT_IMP_;
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}
int s3_deepcopy(void* obj, void** ret){
	try{
		auto surf = static_cast<HM3D::FaceData*>(obj);
		HM3D::FaceData ret_;
		HM3D::DeepCopy(*surf, ret_);
		c2cpp::to_pp(ret_, ret);
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}

int s3_free(void* obj){
	try{
		delete static_cast<HM3D::FaceData*>(obj);
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}
int s3_volume(void* obj, double* ret){
	try{
		auto s = static_cast<HM3D::FaceData*>(obj);
		auto tree = HM3D::Surface::Tree::Assemble(*s);
		HM3D::Surface::R::RevertTree rv(tree);
		*ret = fabs(HM3D::Surface::Volume(*s));
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}
int s3_move(void* obj, double* dx){
	try{
		auto surf = static_cast<HM3D::FaceData*>(obj);
		for (auto p: HM3D::AllVertices(*surf)){
			p->x += dx[0];
			p->y += dx[1];
			p->z += dx[2];
		}
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}
int s3_scale(void* obj, double* pc, double* p0){
	try{
		auto surf = static_cast<HM3D::FaceData*>(obj);
		double xc = pc[0]/100., yc = pc[1]/100., zc = pc[2]/100.;
		for (auto& p: HM3D::AllVertices(*surf)){
			p->x -= p0[0]; p->x *= xc;
			p->y -= p0[1]; p->y *= yc;
			p->z -= p0[2]; p->z *= zc;
		}
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}
int s3_quick_separate(void* obj, int* nret, void*** ret){
	try{
		auto surf = static_cast<HM3D::FaceData*>(obj);
		HM3D::FaceData cp;
		HM3D::DeepCopy(*surf, cp);
		vector<HM3D::FaceData> sep = HM3D::SplitData(cp);
		c2cpp::to_ppp(sep, nret, ret);
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}
int s3_concatenate(int nobjs, void** objs, void** ret){
	try{
		HM3D::FaceData ret_;
		for (auto& it: c2cpp::to_pvec<HM3D::FaceData>(nobjs, objs)){
			ret_.insert(ret_.end(), it->begin(), it->end());
		}
		c2cpp::to_pp(ret_, ret);
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}
int s3_to_hm(void* doc, void* node, void* obj, const char* name, const char* fmt,
		hmcport_callback cb){
	try{
		HMXML::ReaderA* wr = static_cast<HMXML::ReaderA*>(doc);
		HMXML::Reader* sn = static_cast<HMXML::Reader*>(node);
		HM3D::FaceData* c = static_cast<HM3D::FaceData*>(obj);
		HM3D::Export::SurfaceWriter(*c, wr, sn, name, fmt);
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}
