#include "fluent_export_grid3d.hpp"
#include <unordered_map>
#include <fstream>
#include <sstream>
#include "surface_grid3d.hpp"
#include "debug_grid3d.hpp"
#include "serialize_grid3d.hpp"

using namespace HMGrid3D;
namespace hme = HMGrid3D::Export;

HMCallback::FunctionWithCallback<hme::TGridMSH> hme::GridMSH;

namespace{

std::string to_hex(int a){
	std::stringstream os;
	os<<std::hex<<a;
	return os.str();
}

char msh_face_type(HMGrid3D::Face& f){
	int ne = f.n_edges();
	if (ne == 3) return '3';
	if (ne == 4) return '4';
	return '5';
}

char msh_cell_type(HMGrid3D::Cell& c){
	auto dim = c.n_fev();
	int nf = std::get<0>(dim);
	int ne = std::get<1>(dim);
	int nv = std::get<2>(dim);
	//tetrahedral
	if (nf == 4 && ne == 6 && nv == 4)
		return '2';
	//hexahedral
	if (nf == 6 && ne == 12 && nv == 8)
		return '4';
	//pyramid
	if (nf == 5 && ne == 8 && nv == 5)
		return '5';
	//wedge
	if (nf == 5 && ne == 9 && nv == 6)
		return '6';
	//mixed
	return '7';
}

std::map<int, ShpVector<Face>> faces_by_btype(const HMGrid3D::Grid& g){
	std::map<int, ShpVector<Face>> ret;
	//interior zone has always be there
	ret.emplace(std::numeric_limits<int>::min(), ShpVector<Face>());
	for (auto f: g.allfaces()){
		int tp = (f->is_boundary()) ? f->boundary_type : std::numeric_limits<int>::min();
		auto emp = ret.emplace(tp, ShpVector<Face>());
		emp.first->second.push_back(f);
	}
	return ret;
}
std::vector<std::pair<int, int>> left_right_cells(const ShpVector<Face>& data, const ShpVector<HMGrid3D::Cell>& cells){
	std::vector<std::pair<int, int>> ret; ret.reserve(data.size());
	auto _indexer = aa::ptr_container_indexer(cells);
	_indexer.convert();
	for (auto& f: data){
		int p1 = (f->right) ?  _indexer.index(f->right) : -1;
		int p2 = (f->left)  ?  _indexer.index(f->left)  : -1;
		ret.push_back(std::make_pair(p1, p2));
	}

	return ret;
}

std::vector<std::vector<int>> int_face_vertices(const ShpVector<Face>& data, const ShpVector<Vertex>& vert){
	std::vector<std::vector<int>> ret; ret.reserve(data.size());
	auto _indexer = aa::ptr_container_indexer(vert);
	_indexer.convert();
	for (auto& f: data){
		auto sv = f->sorted_vertices();
		ret.push_back(vector<int>());
		ret.back().reserve(sv.size());
		for (auto p: sv){
			ret.back().push_back(_indexer.index(p));
		}
	}
	return ret;
}
char get_common_type(vector<char>::iterator istart, vector<char>::iterator iend){
	for (auto it=istart; it!=iend; ++it) if (*it != *istart) return '0';
	return *istart;
}
typedef std::map<
	std::pair<int, int>,
	std::vector<std::pair<int, int>>
> PeriodicMap;
PeriodicMap assemble_periodic(const std::map<int, ShpVector<Face>>& fzones, const std::map<Face*, Face*>& pfaces){
	std::map<std::pair<int, int>,
		std::vector<std::pair<int, int>>> ret;
	auto perface_index = [&](Face* f){
		int ret = 0;
		auto fnd1 = fzones.find(f->boundary_type);
		for (auto it=fzones.begin(); it!= fnd1; ++it){
			ret += it->second.size();
		}
		auto fnd3 = aa::shp_find(fnd1->second.begin(), fnd1->second.end(), f);
		return ret + fnd3 - (fnd1->second.begin());
	};
	std::map<int, int> btype_zone;
	int zt = 4;
	for (auto it=++fzones.begin(); it!=fzones.end(); ++it){
		btype_zone[it->first] = zt++;
	}
	for (auto& ff: pfaces){
		int perzone = btype_zone[ff.first->boundary_type];
		int shazone = btype_zone[ff.second->boundary_type];
		int perindex = perface_index(ff.first);
		int shaindex = perface_index(ff.second);
		auto emp = ret.emplace(
				std::make_pair(perzone, shazone),
				std::vector<std::pair<int, int>>());
		emp.first->second.emplace_back(perindex, shaindex);
	}
	//sort to always get same result
	for (auto& v: ret) std::sort(v.second.begin(), v.second.end());
	return ret;
}

void zones_names(const PeriodicMap& periodic,
		const std::map<int, ShpVector<Face>>& fzones,
		std::map<int, std::string>& btypes,
		std::map<int, std::string>& bnames){

	std::set<int> periodic_zones, pershadow_zones;
	for (auto& d: periodic){
		int zt1 = d.first.first, zt2 = d.first.second;
		periodic_zones.insert(zt1);
		pershadow_zones.insert(zt2);
	}

	btypes[3] = to_hex(2);
	bnames[3] = "interior";
	int zt = 4;
	for (auto it = ++fzones.begin(); it!=fzones.end(); ++it){
		if (periodic_zones.find(zt) != periodic_zones.end()){
			btypes[zt] = to_hex(12);
			bnames[zt] = "periodic";
		} else if (pershadow_zones.find(zt) != pershadow_zones.end()){
			btypes[zt] = to_hex(8);
			bnames[zt] = "periodic-shadow";
		} else {
			btypes[zt] = to_hex(3);
			bnames[zt] = "wall";
		}
		++zt;
	}
}

/*
struct FaceData{
	int n;                       //total amount of faces
	int istart, iend;            //indicies in total faces array
	int zone_index;              //3,4,5,...
	int zone_type;               //2: interior, ...
	std::string zone_name;       //interior, shadow
	int boundary_index;          //user boundary type
	std::string boundary_name;   //user boundary name
	char face_common_type;       //common type or 0

	static vector<FaceData> to(const SimpleSerialize& ser){
		vector<FaceData> ret;
		_THROW_NOT_IMP_;
		return ret;
	}
	static int zone_by_bnd(const vector<FaceData>& fd, int bnd){
		_THROW_NOT_IMP_;
	}
};

struct PeriodicIntData{
	int zone_index1, zone_index2;
	vector<std::pair<int, int>> periodic;
	struct _dt{
		int b1, b2;
		const Face* f1;
		const Face* f2;
		PeriodicIntData* vptr;
	};

	static vector<PeriodicIntData> to(const SimpleSerialize& ser, const HMGrid3D::Grid& g,
			const vector<FaceData>& face_data,
			const std::map<Face*, Face*>& pface){
		vector<PeriodicIntData> ret; ret.reserve(face_data.size());
		std::vector<_dt> inpvec; inpvec.reserve(pface.size());
		std::map<int, PeriodicIntData*> periodic_bc;
		for (auto& it: pface){
			int b1 = it.first->boundary_type, b2 = it.second->boundary_type;
			auto emp = periodic_bc.emplace(b1, nullptr);
			if (emp.second){
				ret.push_back(PeriodicIntData());
				ret.back().zone_index1 = FaceData::zone_by_bnd(face_data, b1);
				ret.back().zone_index2 = FaceData::zone_by_bnd(face_data, b2);
				emp.first->second = &ret.back();
			}
			inpvec.push_back(PeriodicIntData::_dt {
					b1, b2,
					it.first, it.second,
					emp.first->second});
		}

		//convert Face* to face index
		auto _findexer = aa::shp_container_indexer(std::get<2>(ser._alldata), 
				[](const Face& fc){ return fc.is_boundary(); });
		_findexer.convert();
		for (auto& it: inpvec){
			int find1 = _findexer.index(it.f1);
			int find2 = _findexer.index(it.f2);
			it.vptr->periodic.emplace_back(find1, find2);
		}
		//sort and return
		for (auto& it: ret) std::sort(it.periodic.begin(), it.periodic.end());
		return ret;
	}
};

void gridmsh2(const HMGrid3D::Grid& g, const char* fn,
		std::function<std::string(int)> btype_name,
		std::map<Face*, Face*> pfaces=std::map<Face*, Face*>()){
	if (g.n_cells() == 0) throw std::runtime_error("Exporting blank grid");
	auto& callback = HMCallback::Singleton2::get();

	//Zones:
	//1    - verticies default
	//2    - fluid for cells
	//3    - default interior
	//4..N - bc's faces

	//serialization
	auto ser = SimpleSerialize::Convert(g);
	//... regroup faces ...
	//needed data
	vector<int> face_vert_cells;
	vector<char> cell_types;
	vector<char> face_types;
	char cells_common_type;
	vector<FaceData> face_data = FaceData::to(ser);
	vector<PeriodicIntData> per_data = PeriodicIntData::to(ser, g, face_data, pfaces);

	//=========== Write to file
	callback.silent_step_after(40, "Writing File", 30);
	std::ofstream fs(fn);
	fs.precision(17);
	//header
	fs<<"(0 \"HybMesh to Fluent File\")\n(2 3)\n";

	//Vertices: Zone 1
	callback.subprocess_step_after(10);
	fs<<"(10 (0 1 "<<to_hex(ser.n_vert)<<" 0 3))\n";
	fs<<"(10 (1 1 "<<to_hex(ser.n_vert)<<" 1 3)(\n";
	for (int i=0; i<ser.vert.size(); i+=3){
		fs<<ser.vert[i]<<" "<<ser.vert[i+1]<<" "<<ser.vert[i+2]<<"\n";
	}
	fs<<"))\n";

	//Cells: Zone 2
	callback.subprocess_step_after(10);
	fs<<"(12 (0 1 "<<to_hex(ser.n_cells)<<" 0))\n";
	fs<<"(12 (2 1 "<<to_hex(ser.n_cells)<<" 1 "<<cells_common_type<<")";
	if (cells_common_type != '0') fs<<")\n";
	else {
		fs<<"(\n";
		for (auto s: cell_types) fs<<s<<" ";
		fs<<"\n))\n";
	}

	//Faces: Zones 3+it
	callback.subprocess_step_after(10);
	fs<<"(13 (0 1 "<<to_hex(ser.n_faces)<<" 0))\n";

	int t=0;
	for (int k=0; k<face_data.size(); ++k){
		auto& fd = face_data[k];
		if (fd.n == 0) continue;
		fs<<"(13 ("<<to_hex(fd.zone_index)<<" "<<to_hex(fd.istart+1)<<" "<<to_hex(fd.iend)<<" ";
		fs<<to_hex(fd.zone_type)<<" "<<fd.face_common_type<<")(\n";
		bool write_sz = (fd.face_common_type == '0' || fd.face_common_type == '5');

		for (int i=fd.istart; i<fd.iend; ++i){
			int sz = face_vert_cells[t++];
			if (write_sz) fs<<to_hex(sz)<<" ";
			for (int j=0; j<sz; ++j) fs<<to_hex(face_vert_cells[t++] + 1)<<" ";
			fs<<to_hex(face_vert_cells[t++] + 1)<<" ";
			fs<<to_hex(face_vert_cells[t++] + 1)<<"\n";
		}
	}

	//Periodic
	t = 0;
	for (auto& pd: per_data){
		int sz = pd.periodic.size();
		fs<<"(18 ("<<to_hex(t+1)<<" "<<to_hex(t+sz)<<" ";
		fs<<to_hex(pd.zone_index1)<<" "<<to_hex(pd.zone_index2)<<")(\n";
		for (int i=0; i<sz; ++i){
			fs<<to_hex(pd.periodic[i].first+1)<<" ";
			fs<<to_hex(pd.periodic[i].second+1)<<"\n";
		}
		fs<<"))\n";
		t += sz;
	}

	//Boundary features
	fs<<"(45 (2 fluid fluid)())\n";
	for (auto& fd: face_data){
		fs<<"(45 ("<<fd.zone_index<<" "<<fd.zone_name<<" "<<fd.boundary_name<<")())\n";
	}

	callback.fin();
	
}
*/

//save to fluent main function
void gridmsh(HMCallback::Caller2& callback, const HMGrid3D::Grid& g, std::string fn,
		hme::BFun btype_name,
		std::map<Face*, Face*> pfaces=std::map<Face*, Face*>()){
	if (g.n_cells() == 0) throw std::runtime_error("Exporting blank grid");
	//Zones:
	//1    - verticies default
	//2    - fluid for cells
	//3    - default interior
	//4..N - bc's faces

	// ===== Needed data
	callback.silent_step_after(60, "Assembling connectivity", 70);
	//* vertices table
	callback.subprocess_step_after(10);
	auto vert = g.allvertices();
	//* faces by boundary_type grouped by zones
	callback.subprocess_step_after(10);
	std::map<int, ShpVector<Face>> facezones = faces_by_btype(g);
	//* faces as they would be in resulting file
	callback.subprocess_step_after(10);
	ShpVector<Face> allfaces;
	for (auto& it: facezones) allfaces.insert(allfaces.end(), it.second.begin(), it.second.end());
	//* left/right cell indicies for each face in integer representation
	callback.subprocess_step_after(10);
	std::vector<std::pair<int, int>> face_adjacents = left_right_cells(allfaces, g.allcells()); 
	//* face->vertices connectivity in integer representation
	callback.subprocess_step_after(10);
	std::vector<std::vector<int>> face_vertices = int_face_vertices(allfaces, vert);
	//* cells, faces types for each entry;
	callback.subprocess_step_after(10);
	std::vector<char> ctypes;   //cell types array
	std::vector<char> ftypes;   //face types for each zone
	for (auto& c: g.allcells()) ctypes.push_back(msh_cell_type(*c));
	for (auto& f: allfaces)     ftypes.push_back(msh_face_type(*f));
        //* common cell types or '0' if they differ, common face type for each zone
	char cell_common_type = get_common_type(ctypes.begin(), ctypes.end());
	std::vector<char> facezones_common_type;
	auto ftypes_it = ftypes.begin();
	for (auto& v: facezones){
		auto itnext = ftypes_it + v.second.size();
		facezones_common_type.push_back(get_common_type(ftypes_it, itnext));
		ftypes_it = itnext;
	}
	//* <periodic zonetype, shadow zonetype> -> vector of periodic face, shadow face>
	callback.subprocess_step_after(10);
	PeriodicMap periodic_data = assemble_periodic(facezones, pfaces);
	//* name and type of each zone by zone id: 2-interior, 3-wall, 8-periodic-shadow, 12-periodic
	std::map<int, std::string> facezones_btype, facezones_bname;
	zones_names(periodic_data, facezones, facezones_btype, facezones_bname);
	callback.subprocess_fin();
	
	//=========== Write to file
	callback.silent_step_after(40, "Writing File", 30);
	std::ofstream fs(fn);
	fs.precision(17);
	//header
	fs<<"(0 \"HybMesh to Fluent File\")\n(2 3)\n";

	//Vertices: Zone 1
	callback.subprocess_step_after(10);
	fs<<"(10 (0 1 "<<to_hex(vert.size())<<" 0 3))\n";
	fs<<"(10 (1 1 "<<to_hex(vert.size())<<" 1 3)(\n";
	for (auto p: vert){
		fs<<p->x<<" "<<p->y<<" "<<p->z<<"\n";
	}
	fs<<"))\n";

	//Cells: Zone 2
	callback.subprocess_step_after(10);
	fs<<"(12 (0 1 "<<to_hex(g.n_cells())<<" 0))\n";
	fs<<"(12 (2 1 "<<to_hex(g.n_cells())<<" 1 "<<cell_common_type<<")";
	if (cell_common_type != '0') fs<<")\n";
	else {
		fs<<"(\n";
		for (auto s: ctypes) fs<<s<<" ";
		fs<<"\n))\n";
	}

	//Faces: Zones 3+it
	callback.subprocess_step_after(10);
	fs<<"(13 (0 1 "<<to_hex(allfaces.size())<<" 0))\n";
	int zonetype = 3;
	int iface = 0;
	int it = 0;
	for (auto& fz: facezones){
		if (fz.second.size() == 0) { ++it; ++zonetype; continue;}
		std::string first_index = to_hex(iface+1);
		std::string last_index = to_hex(iface+fz.second.size());
		fs<<"(13 ("<<to_hex(zonetype)<<" "<<first_index<<" "<<last_index<<" ";
		fs<<facezones_btype[zonetype]<<" "<<facezones_common_type[it]<<")(\n";
		for (auto f: fz.second){
			if (facezones_common_type[it] == '0' ||
			    facezones_common_type[it] == '5')
				fs<<to_hex(face_vertices[iface].size())<<" ";
			for (auto i: face_vertices[iface]) fs<<to_hex(i+1)<<" ";
			fs<<to_hex(face_adjacents[iface].first+1)<<" "<<to_hex(face_adjacents[iface].second+1);
			fs<<"\n";
			++iface;
		}
		fs<<"))\n";
		++it; ++zonetype;
	}
	//Periodic
	if (pfaces.size() > 0){
		int it = 0;
		for (auto& pd: periodic_data){
			int sz = pd.second.size();
			fs<<"(18 ("<<to_hex(it+1)<<" "<<to_hex(it+sz)<<" ";
			fs<<to_hex(pd.first.first)<<" "<<to_hex(pd.first.second)<<")(\n";
			for (auto& v: pd.second){
				fs<<to_hex(v.first+1)<<" "<<to_hex(v.second+1)<<"\n";
			}
			fs<<"))\n";
			it += sz;
		}
	}

	//Boundary features
	fs<<"(45 (2 fluid fluid)())\n";
	fs<<"(45 (3 interior default-interior)())\n";
	zonetype = 4;
	for (auto it = std::next(facezones.begin()); it!=facezones.end(); ++it){
		fs<<"(45 ("<<zonetype<<" "<<facezones_bname[zonetype]<<" "<<btype_name(it->first)<<")())\n";
		++zonetype;
	}
}

}

void hme::PeriodicDataEntry::assemble(HMGrid3D::Grid& g, std::map<Face*, Face*>& outmap){
	HMGrid3D::Surface periodic_surf = HMGrid3D::Surface::FromBoundaryType(g, bt, 1);
	int dir = (reversed) ? -1 : 1;
	HMGrid3D::Surface shadow_surf = HMGrid3D::Surface::FromBoundaryType(g, bt_shadow, dir);
	//this (along with not null dir) guaranties that all bondary edges
	//would have same direction that is needed to call ExtractBoundary
	for (auto f: periodic_surf.allfaces()) f->correct_edge_directions();
	for (auto f: shadow_surf.allfaces()) f->correct_edge_directions();

	//boundary
	ShpVector<Edge> bnd1 = Surface::ExtractBoundary(periodic_surf, v);
	ShpVector<Edge> bnd2 = Surface::ExtractBoundary(shadow_surf, v_shadow);
	//boundary sizes should be equal
	if (bnd1.size() != bnd2.size() || bnd1.size() == 0)
		throw std::runtime_error("Periodic merging failed");

	//subsurface bounded by bnd1, bnd2 (if initial surface are multiply connected)
	periodic_surf = HMGrid3D::Surface::SubSurface(periodic_surf, bnd1[0]->vertices[0].get());
	shadow_surf = HMGrid3D::Surface::SubSurface(shadow_surf, bnd2[0]->vertices[0].get());

	//Rearranging to match topology
	HMGrid3D::Surface::FaceRearrange(periodic_surf, bnd1[0].get());
	HMGrid3D::Surface::FaceRearrange(shadow_surf, bnd2[0].get());

	//check topology: if ok then map surface else throw
	if (HMGrid3D::Surface::MatchTopology(periodic_surf, shadow_surf)){
		auto it1 = periodic_surf.faces.begin();
		auto it2 = shadow_surf.faces.begin();
		while (it1 != periodic_surf.faces.end()){
			outmap[(*it1++).get()] = (*it2++).get();
		}
	} else {
		throw std::runtime_error("Periodic merging failed");
	}
}

HMGrid3D::Grid hme::PeriodicData::assemble(const HMGrid3D::Grid& g,
		std::map<Face*, Face*>& outmap){
	if (size() == 0) return g;
	HMGrid3D::Grid ret = HMGrid3D::Constructor::Copy::ShallowVertices(g);
	for (auto& d: data) d.assemble(ret, outmap);
	return ret;
}

void hme::TGridMSH::_run(const Grid& g, std::string fn,
		BFun btype_name, PeriodicData periodic){
	callback.step_after(30, "Periodic merging");
	std::map<Face*, Face*> periodic_cells;
	HMGrid3D::Grid gp = periodic.assemble(g, periodic_cells);
	return gridmsh(callback, gp, fn, btype_name, periodic_cells);
}

void hme::TGridMSH::_run(const HMGrid3D::Grid& g, std::string fn,
		PeriodicData periodic){
	return _run(g, fn, def_bfun, periodic);
}

void hme::TGridMSH::_run(const HMGrid3D::Grid& g, std::string fn, BFun btype_name){
	return gridmsh(callback, g, fn, btype_name);
}

void hme::TGridMSH::_run(const HMGrid3D::Grid& g, std::string fn){
	return _run(g, fn, def_bfun);
}


