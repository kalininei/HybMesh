#include "fluent_export3d.hpp"
#include <unordered_map>
#include <fstream>
#include <sstream>
#include "surface.hpp"
#include "debug3d.hpp"
#include "serialize3d.hpp"

using namespace HM3D;
namespace hme = HM3D::Export;

HMCallback::FunctionWithCallback<hme::TGridMSH> hme::GridMSH;

namespace{

std::string to_hex(int a){
	std::stringstream os;
	os<<std::hex<<a;
	return os.str();
}

char msh_face_type(HM3D::Face& f){
	int ne = f.edges.size();
	if (ne == 3) return '3';
	if (ne == 4) return '4';
	return '5';
}

char msh_cell_type(HM3D::Cell& c){
	auto nfev = c.n_fev();
	int nf = std::get<0>(nfev);
	int ne = std::get<1>(nfev);
	int nv = std::get<2>(nfev);
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

std::map<int, ShpVector<Face>> faces_by_btype(const GridData& g){
	std::map<int, ShpVector<Face>> ret;
	//interior zone has always be there
	ret.emplace(std::numeric_limits<int>::min(), ShpVector<Face>());
	for (auto f: g.vfaces){
		int tp = (f->is_boundary()) ? f->boundary_type : std::numeric_limits<int>::min();
		auto emp = ret.emplace(tp, ShpVector<Face>());
		emp.first->second.push_back(f);
	}
	return ret;
}
std::vector<std::pair<int, int>> left_right_cells(const ShpVector<Face>& data, const ShpVector<HM3D::Cell>& cells){
	std::vector<std::pair<int, int>> ret; ret.reserve(data.size());
	auto _indexer = aa::ptr_container_indexer(cells);
	_indexer.convert();
	for (auto& f: data){
		int p1 = (f->has_right_cell()) ?  _indexer.index(f->right.lock()) : -1;
		int p2 = (f->has_left_cell())  ?  _indexer.index(f->left.lock())  : -1;
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

//save to fluent main function
void gridmsh(HMCallback::Caller2& callback, const GridData& g, std::string fn,
		hme::BFun btype_name,
		std::map<Face*, Face*> pfaces=std::map<Face*, Face*>()){
	if (g.vcells.size() == 0) throw std::runtime_error("Exporting blank grid");
	//Zones:
	//1    - verticies default
	//2    - fluid for cells
	//3    - default interior
	//4..N - bc's faces

	// ===== Needed data
	callback.silent_step_after(60, "Assembling connectivity", 70);
	//* vertices table
	callback.subprocess_step_after(10);
	auto vert = g.vvert;
	//* faces by boundary_type grouped by zones
	callback.subprocess_step_after(10);
	std::map<int, ShpVector<Face>> facezones = faces_by_btype(g);
	//* faces as they would be in resulting file
	callback.subprocess_step_after(10);
	ShpVector<Face> allfaces;
	for (auto& it: facezones) allfaces.insert(allfaces.end(), it.second.begin(), it.second.end());
	//* left/right cell indicies for each face in integer representation
	callback.subprocess_step_after(10);
	std::vector<std::pair<int, int>> face_adjacents = left_right_cells(allfaces, g.vcells); 
	//* face->vertices connectivity in integer representation
	callback.subprocess_step_after(10);
	std::vector<std::vector<int>> face_vertices = int_face_vertices(allfaces, vert);
	//* cells, faces types for each entry;
	callback.subprocess_step_after(10);
	std::vector<char> ctypes;   //cell types array
	std::vector<char> ftypes;   //face types for each zone
	for (auto& c: g.vcells) ctypes.push_back(msh_cell_type(*c));
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
	fs<<"(12 (0 1 "<<to_hex(g.vcells.size())<<" 0))\n";
	fs<<"(12 (2 1 "<<to_hex(g.vcells.size())<<" 1 "<<cell_common_type<<")";
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
	fs<<"(45 (2 fluid fluid 1)())\n";
	fs<<"(45 (3 interior default-interior 1)())\n";
	zonetype = 4;
	for (auto it = std::next(facezones.begin()); it!=facezones.end(); ++it){
		fs<<"(45 ("<<zonetype<<" "<<facezones_bname[zonetype]<<" "<<btype_name(it->first)<<" 1)())\n";
		++zonetype;
	}
}

}

void hme::PeriodicDataEntry::assemble(GridData& g, std::map<Face*, Face*>& outmap){
	auto surfs = HM3D::Surface::GridSurfaceBType(g);
	HM3D::FaceData periodic_surf = surfs[bt];
	HM3D::FaceData shadow_surf = surfs[bt_shadow];
	
	bool dir2 = (reversed) ? true : false;
	Surface::RevertGridSurface r1(periodic_surf, false);
	Surface::RevertGridSurface r2(shadow_surf, dir2);
	r1.make_permanent();
	r2.make_permanent();

	//this (along with not null dir) guaranties that all bondary edges
	//would have same direction that is needed to call ExtractBoundary
	for (auto f: periodic_surf) f->correct_edge_directions();
	for (auto f: shadow_surf) f->correct_edge_directions();

	//boundary
	ShpVector<Edge> bnd1 = Surface::ExtractBoundary(periodic_surf, v);
	ShpVector<Edge> bnd2 = Surface::ExtractBoundary(shadow_surf, v_shadow);
	//boundary sizes should be equal
	if (bnd1.size() != bnd2.size() || bnd1.size() == 0)
		throw std::runtime_error("Periodic merging failed");

	//subsurface bounded by bnd1, bnd2 (if initial surface are multiply connected)
	periodic_surf = HM3D::Surface::SubSurface(periodic_surf, bnd1[0]->vertices[0].get());
	shadow_surf = HM3D::Surface::SubSurface(shadow_surf, bnd2[0]->vertices[0].get());

	//Rearranging to match topology
	HM3D::Surface::FaceRearrange(periodic_surf, bnd1[0].get());
	HM3D::Surface::FaceRearrange(shadow_surf, bnd2[0].get());

	//check topology: if ok then map surface else throw
	if (HM3D::Surface::MatchTopology(periodic_surf, shadow_surf)){
		auto it1 = periodic_surf.begin();
		auto it2 = shadow_surf.begin();
		while (it1 != periodic_surf.end()){
			outmap[(*it1++).get()] = (*it2++).get();
		}
	} else {
		throw std::runtime_error("Periodic merging failed");
	}
}

HM3D::GridData hme::PeriodicData::assemble(const GridData& g,
		std::map<Face*, Face*>& outmap){
	if (size() == 0) return g;
	HM3D::GridData ret;
	HM3D::DeepCopy(g, ret, 2);
	for (auto& d: data) d.assemble(ret, outmap);
	return ret;
}

void hme::TGridMSH::_run(const Ser::Grid& g, std::string fn,
		BFun btype_name, PeriodicData periodic){
	callback->step_after(30, "Periodic merging");
	std::map<Face*, Face*> periodic_cells;
	GridData gp = periodic.assemble(g.grid, periodic_cells);
	return gridmsh(*callback, gp, fn, btype_name, periodic_cells);
}

void hme::TGridMSH::_run(const HM3D::Ser::Grid& g, std::string fn,
		PeriodicData periodic){
	return _run(g, fn, def_bfun, periodic);
}

void hme::TGridMSH::_run(const HM3D::Ser::Grid& g, std::string fn, BFun btype_name){
	return gridmsh(*callback, g.grid, fn, btype_name);
}

void hme::TGridMSH::_run(const HM3D::Ser::Grid& g, std::string fn){
	return _run(g, fn, def_bfun);
}
