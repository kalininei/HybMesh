#include "export3d_gmsh.hpp"
#include <sstream>
#include <fstream>
#include "addalgo.hpp"
#include "serialize3d.hpp"
#include "export3d_vtk.hpp"
#include "surface.hpp"

using namespace HM3D;
namespace hme = HM3D::Export;

HMCallback::FunctionWithCallback<hme::TGridGMSH> hme::GridGMSH;

void hme::TGridGMSH::_run(const Ser::Grid& ser, std::string fn, BFun bfun){
	const GridData& grid = ser.grid;
	callback->step_after(25, "Faces assembling");
	auto fv = ser.face_vertex();
	//face data
	std::map<int, FaceData> srfs = Surface::GridSurfaceBType(grid);
	std::map<int, vector<vector<int>>> psrfs;
	int totfaces = 0;
	aa::enumerate_ids_pvec(grid.vfaces);
	for (auto& s: srfs){
		auto er = psrfs.emplace(s.first, vector<vector<int>>(s.second.size()));
		auto& vv = er.first->second;
		totfaces+=s.second.size();
		//enumerate nodes so that all cells be on their left side
		for (int i=0; i<s.second.size(); ++i){
			auto fc = s.second[i];
			int iface = fc->id;
			vv[i] = fv[iface];
			if (!fc->has_left_cell()) std::reverse(vv[i].begin()+1, vv[i].end());
		}
	}
	callback->step_after(30, "Cells assembling");
	//cell data
	auto cp = hme::vtkcell_expression::cell_assembler(ser, fv);
	int ientity = psrfs.rbegin()->first+1;
	auto cell_line = [&cp, &ientity](int num)->std::string{
		auto& c = cp[num];
		std::ostringstream os;
		int tp;
		switch (c.celltype){
			case 10: tp = 4; break;
			case 12: tp = 5; break;
			case 14: tp = 7; break;
			case 13: tp = 6; break;
			default: assert(false);
		}
		os<<num+1<<" "<<tp<<" 2 "<<ientity<<" "<<ientity;
		if (tp != 6){
			for (size_t i=0; i<c.pts.size(); ++i)
				os<<" "<<c.pts[i]+1;
		} else {
			os<<" "<<c.pts[0]+1<<" "<<c.pts[2]+1<<" "<<c.pts[1]+1;
			os<<" "<<c.pts[3]+1<<" "<<c.pts[5]+1<<" "<<c.pts[4]+1;
		}
		return os.str();
	};

	std::ofstream of(fn);
	of.precision(17);
	//header
	of<<"$MeshFormat"<<std::endl;
	of<<"2.2 0 8"<<std::endl;
	of<<"$EndMeshFormat"<<std::endl;
	//physical tags
	of<<"$PhysicalNames"<<std::endl;
	of<<psrfs.size()+1<<std::endl;
	for(auto& s: psrfs){
		of<<"2 "<<s.first<<" \""<<bfun(s.first)<<"\""<<std::endl;
	}
	of<<"3 "<<ientity<<" \"interior\""<<std::endl;
	of<<"$EndPhysicalNames"<<std::endl;
	//nodes
	callback->step_after(25, "Nodes writing");
	of<<"$Nodes"<<std::endl;
	{
		of<<ser.n_vert()<<std::endl;
		auto& vert = ser.vert();
		for (int i=0; i<ser.n_vert(); ++i){
			of<<i+1<<" "<<vert[3*i]<<" "<<vert[3*i+1]<<" "<<vert[3*i+2]<<std::endl;
		}
	}
	of<<"$EndNodes"<<std::endl;
	of<<"$Elements"<<std::endl;
	of<<ser.n_cells()+totfaces<<std::endl;
	//cells
	{
		callback->step_after(10, "Interior cells");
		//3d cells
		for (int i=0; i<ser.n_cells(); ++i){
			of<<cell_line(i)<<std::endl;
		}
		int icell=ser.n_cells();
		callback->step_after(10, "Boundary cells");
		//boundary
		for (auto& s: psrfs){
			int entity = s.first;
			for (auto& k: s.second){
				int tp = (k.size()==3)?2:3;
				of<<++icell<<" "<<tp<<" 2 "<<entity<<" "<<entity;
				for (auto& n: k) of<<" "<<n+1;
				of<<std::endl;
			}
		}
	}
	of<<"$EndElements"<<std::endl;
}

void hme::TGridGMSH::_run(const Ser::Grid& ser, std::string fn){
	return _run(ser, fn, def_bfun);
}
