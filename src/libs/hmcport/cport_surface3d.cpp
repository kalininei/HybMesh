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
int s3_tab_vertices(void* obj, double* ret){
	try{
		auto surf = static_cast<HM3D::FaceData*>(obj);
		for (auto& v: HM3D::AllVertices(*surf)){
			*ret++ = v->x;
			*ret++ = v->y;
			*ret++ = v->z;
		}
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}
int s3_tab_edgevert(void* obj, int* ret){
	try{
		auto surf = static_cast<HM3D::FaceData*>(obj);
		auto ae = AllEdges(*surf);
		auto av = AllVertices(*surf);
		aa::enumerate_ids_pvec(av);
		for (auto& e: ae){
			*ret++ = e->pfirst()->id;
			*ret++ = e->plast()->id;
		}
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}
int s3_tab_facedim(void* obj, int* ret){
	try{
		auto surf = static_cast<HM3D::FaceData*>(obj);
		for (auto& f: *surf){
			*ret++ = f->edges.size();
		}
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}
int s3_tab_faceedge(void* obj, int* nret, int** cret){
	try{
		auto surf = static_cast<HM3D::FaceData*>(obj);
		*nret = 0;
		for (auto& f: *surf) *nret+=f->edges.size();
		*cret = new int[*nret];
		int* ret = *cret;
		auto ae = HM3D::AllEdges(*surf);
		aa::enumerate_ids_pvec(ae);
		for (auto& f: *surf){
			for (auto& e: f->edges){
				*ret++ = e->id;
			}
		}
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}
int s3_tab_facevert(void* obj, int* nret, int** cret){
	try{
		auto surf = static_cast<HM3D::FaceData*>(obj);
		*nret = 0;
		for (auto& f: *surf) *nret+=f->edges.size();
		*cret = new int[*nret];
		int* ret = *cret;
		auto av = HM3D::AllVertices(*surf);
		aa::enumerate_ids_pvec(av);
		for (auto& f: *surf){
			for(auto& v: f->sorted_vertices()){
				*ret++ = v->id;
			}
		}
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}
int s3_tab_centers(void* obj, double* ret){
	try{
		auto surf = static_cast<HM3D::FaceData*>(obj);
		for (auto& f: *surf){
			Point3 pnt(0, 0, 0);
			for (auto& v: AllVertices(f->edges)){
				pnt += *v;
			}
			pnt /= f->edges.size();
			*ret++ = pnt.x;
			*ret++ = pnt.y;
			*ret++ = pnt.z;
		}
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

int s3_assign_boundary_types(void* obj, int* bnd, int** revdif){
	try{
		auto surf = static_cast<HM3D::FaceData*>(obj);
		auto whole_assign = [&surf](int bt, std::map<int, int>& mp){
			for (int k=0; k<surf->size(); ++k){
				mp[k] = bt;
			}
		};
		auto bt_by_index = [&surf](int ind)->int&{
			return (*surf)[ind]->boundary_type;
		};
		c2cpp::assign_boundary_types(bnd, revdif, whole_assign, bt_by_index);
		return HMSUCCESS;
	} catch (std::exception& e){
		add_error_message(e.what());
		return HMERROR;
	}
}
