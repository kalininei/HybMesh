#include "export2d_fluent.hpp"
#include <sstream>
#include <fstream>
#include <unordered_map>
#include "algos.hpp"
#include "cont_assembler.hpp"

using namespace HM2D;
namespace hme=HM2D::Export;

namespace{

std::string default_bfun(int i){
	return std::string("boundary") + std::to_string(i);
}

std::string to_hex(int a){
	std::stringstream os;
	os<<std::hex<<a;
	return os.str();
}

struct FaceData{
	int n;                       //total amount of faces
	int istart, iend;            //indicies in total faces array
	int zone_index;              //3,4,5,...
	std::string zone_name;       //interior, shadow
	int zone_type;               //2: interior, 3: wall, 12/8: periodic/shadow
	int boundary_index;          //user boundary type
	std::string boundary_name;   //user boundary name

	static FaceData* zone_by_bindex(vector<FaceData>& fd, int bindex){
		for (auto& v: fd) { if (v.boundary_index == bindex) return &v; }
		assert(false);
		return 0;
	}

	static vector<FaceData> to(GridData& gg, hme::BNamesFun bnames, hme::PeriodicData pd){
		// ==== renumber edges according to boundary types
		std::map<int, vector<int>> bmap;
		for (int i=0; i<gg.vedges.size(); ++i){
			int bt = std::numeric_limits<int>::min();
			if (gg.vedges[i]->is_boundary()) bt = gg.vedges[i]->boundary_type;
			auto emp = bmap.emplace(bt, vector<int>());
			emp.first->second.push_back(i);
		}
		EdgeData new_edges; new_edges.reserve(new_edges.size());
		for (auto& m: bmap){
			for (int ind: m.second){
				new_edges.push_back(gg.vedges[ind]);
			}
		}
		std::swap(gg.vedges, new_edges);

		// ==== fill return array
		vector<FaceData> ret(1);
		int icur = 0;
		//internal faces
		FaceData* fdcur = &ret.back();
		fdcur->istart = 0;
		fdcur->zone_index = 3;
		fdcur->zone_name = "interior";
		fdcur->zone_type = 2;
		fdcur->boundary_name = "default-interior";
		while (!gg.vedges[icur]->is_boundary() && icur<gg.vedges.size()) ++icur;
		fdcur->iend = icur;
		fdcur->n = fdcur->iend - fdcur->istart;
		//boundary faces
		while (icur < gg.vedges.size()){
			ret.emplace_back();
			fdcur = &ret.back();
			fdcur->istart = icur;
			fdcur->boundary_index = gg.vedges[icur]->boundary_type;
			fdcur->boundary_name = bnames(fdcur->boundary_index);
			bool fnd_per = std::find(pd.b1.begin(), pd.b1.end(), fdcur->boundary_index) != pd.b1.end();
			bool fnd_shd = std::find(pd.b2.begin(), pd.b2.end(), fdcur->boundary_index) != pd.b2.end();
			if (fnd_per){
				fdcur->zone_name = "periodic";
				fdcur->zone_type = 12;
			}else if (fnd_shd){
				fdcur->zone_name = "periodic-shadow";
				fdcur->zone_type = 8;
			} else {
				fdcur->zone_name = "wall";
				fdcur->zone_type = 3;
			}
			fdcur->zone_index = (++ret.rbegin())->zone_index + 1;
			while(icur<gg.vedges.size() && gg.vedges[icur]->boundary_type == fdcur->boundary_index) ++icur;
			fdcur->iend = icur;
			fdcur->n = fdcur->iend - fdcur->istart;
		}
		return ret;
	}
};
char get_common_type(vector<char>::iterator istart, vector<char>::iterator iend){
	for (auto it=istart; it!=iend; ++it) if (*it != *istart) return '0';
	return *istart;
}
char msh_cell_type(const Cell& c){
	int ne = c.edges.size();
	if (ne == 3) return '1';
	if (ne == 4) return '3';
	throw std::runtime_error("invalid cell type for Fluent export");
}

}

std::vector<int> hme::PeriodicData::assemble(GridData& g) const {
	if (size() == 0) return {};
	//all periodic boundary types to single array
	std::vector<int> btypes;
	btypes.insert(btypes.end(), b1.begin(), b1.end());
	btypes.insert(btypes.end(), b2.begin(), b2.end());
	//sort out only those faces which are boundary
	std::vector<Edge*> bnd_edges;
	std::vector<int> bnd_edges_indicies;
	std::vector<int> bnd_edges_btp;
	for (auto i=0; i<g.vedges.size(); ++i) if (g.vedges[i]->is_boundary()){
		int bt = g.vedges[i]->boundary_type;
		if (std::find(btypes.begin(), btypes.end(), bt) != btypes.end()){
			bnd_edges.push_back(g.vedges[i].get());
			bnd_edges_indicies.push_back(i);
			bnd_edges_btp.push_back(bt);
		}
	}
	//revert edges. After this == for grid edges will not be reliable.
	//However we are not goint to use it further
	for (int i=0; i<bnd_edges.size(); ++i){
		int bt = bnd_edges_btp[i];
		Edge* e = bnd_edges[i];
		int index_in_pd;
		bool is_per;
		auto fnd1 = std::find(b1.begin(), b1.end(), bt);
		if (fnd1 != b1.end()){
			index_in_pd = fnd1 - b1.begin();
			is_per = true;
		} else {
			auto fnd2 = std::find(b2.begin(), b2.end(), bt);
			assert(fnd2 != b2.end());
			index_in_pd = fnd2 - b2.begin();
			is_per = false;
		}
		//periodic faces have cells on their left
		//shadow faces have cells on their right/left for rev = true/false
		bool need_rev = false;
		if (is_per && e->no_left_cell()) need_rev = true;
		else if (isrev[index_in_pd]){
			if (e->no_right_cell()) need_rev = true;
		} else {
			if (e->no_left_cell()) need_rev = true;
		}
		if (need_rev){
			e->reverse();
		}
	}
	//build contour collection from boundary grid edges
	//all edges are directed so that cells are on left side
	std::vector<HM2D::EdgeData> ecols;
	for (int i=0; i<btypes.size(); ++i){
		ecols.emplace_back();
		for (int k=0; k<bnd_edges.size(); ++k) if (bnd_edges_btp[k] == btypes[i]){
			Edge* ge = bnd_edges[k];
			shared_ptr<HM2D::Edge> ce( new HM2D::Edge(0, 0));
			ce->vertices[0] = std::make_shared<HM2D::Vertex>(*ge->first());
			ce->vertices[1] = std::make_shared<HM2D::Vertex>(*ge->last());
			if (ge->no_left_cell()) ce->reverse();
			ecols.back().push_back(ce);
		}
		HM2D::ECol::Algos::MergePoints(ecols.back());
	}
	//assemble contours
	std::vector<HM2D::EdgeData> assembled_b_contours;
	for (int i=0; i<ecols.size(); ++i){
		auto ar = HM2D::Contour::Assembler::AllContours(ecols[i]);
		if (ar.size() != 1){
			throw std::runtime_error("Periodic merging failed");
		}
		//check direction
		if (ar[0].size() > 1 && ar[0][0]->last() != ar[0][1]->first()) HM2D::Contour::Reverse(ar[0]);
		assembled_b_contours.push_back(std::move(ar[0]));
	}
	//map for adressing
	std::unordered_map<HM2D::Edge*, int> e_cont_to_grid;
	for (auto& it1: assembled_b_contours)
	for (auto& it2: it1){
		//int p1 = static_cast<GridPoint*>(it2->pstart)->get_ind();
		//int p2 = static_cast<GridPoint*>(it2->pend)->get_ind();
		auto v1 = it2->first().get();
		auto v2 = it2->last().get();
		int gi;
		for (gi=0; gi<bnd_edges.size(); ++gi){
			auto e = bnd_edges[gi];
			if (*e->first() != *v1 && *e->first() != *v2) continue;
			if (*e->last() != *v1 && *e->last() != *v2) continue;
			break;
		}
		assert(gi<bnd_edges.size());
		e_cont_to_grid.emplace(it2.get(), bnd_edges_indicies[gi]);
	}
	//build relations
	std::vector<HM2D::Edge*> cont_edges_relations;
	for (int k=0; k<size(); ++k){
		int k1 = std::find(btypes.begin(), btypes.end(), b1[k]) - btypes.begin();
		int k2 = std::find(btypes.begin(), btypes.end(), b2[k]) - btypes.begin();
		auto c1 = &assembled_b_contours[k1];
		auto c2 = &assembled_b_contours[k2];
		if (c1->size() != c2->size() ||
				HM2D::Contour::IsClosed(*c1) ||
				HM2D::Contour::IsClosed(*c2)){
			throw std::runtime_error("Periodic merging failed");
		}
		if (isrev[k]) HM2D::Contour::Reverse(*c2);
		for (int i=0; i<c1->size(); ++i){
			cont_edges_relations.push_back((*c1)[i].get());
			cont_edges_relations.push_back((*c2)[i].get());
		}
	}
	// assemble result
	std::vector<int> ret; ret.reserve(cont_edges_relations.size());
	for (auto e: cont_edges_relations) ret.push_back(e_cont_to_grid[e]);
	return ret;
}

void hme::GridMSH(const GridData& gg, std::string fn, hme::BNamesFun bnames, PeriodicData pd){
	//Zones:
	//1    - verticies default
	//2    - fluid for cells
	//3    - default interior
	//4..N - bc's faces
	
	//make a deep copy of grid since some data order changes will be made.
	GridData g;
	HM2D::DeepCopy(gg, g);

	//Needed data
	vector<char> ctypes(g.vcells.size());
	for (int i=0; i<g.vcells.size(); ++i) ctypes[i] = msh_cell_type(*g.vcells[i]);
	char cell_common_type = get_common_type(ctypes.begin(), ctypes.end());

	//edges will be resorted
	vector<FaceData> face_data;
	face_data = FaceData::to(g, bnames, pd);
	//edges will be reverted if necessary
	vector<int> periodic_relations = pd.assemble(g);

	std::ofstream fs(fn);
	fs.precision(16);
	//header
	fs<<"(0 \"HybMesh to Fluent File\")\n(2 2)\n";

	//Vertices: Zone 1
	fs<<"(10 (0 1 "<<to_hex(g.vvert.size())<<" 0 2))\n";
	fs<<"(10 (1 1 "<<to_hex(g.vvert.size())<<" 1 2)(\n";
	for (int i=0; i<g.vvert.size(); ++i){
		auto& p = g.vvert[i];
		fs<<p->x<<" "<<p->y<<"\n";
	}
	fs<<"))\n";

	//Cells: Zone2
	fs<<"(12 (0 1 "<<to_hex(g.vcells.size())<<" 0))\n";
	fs<<"(12 (2 1 "<<to_hex(g.vcells.size())<<" 1 "<<cell_common_type<<")";
	if (cell_common_type != '0') fs<<")\n";
	else {
		fs<<"(\n";
		for (auto s: ctypes) fs<<s<<" ";
		fs<<"\n))\n";
	}

	//Faces: Zones 3+it
	fs<<"(13 (0 1 "<<to_hex(g.vedges.size())<<" 0))\n";
	aa::enumerate_ids_pvec(g.vcells);
	aa::enumerate_ids_pvec(g.vvert);
	int iface = 0;
	for (auto& fd: face_data){
		if (fd.n == 0) continue;
		fs<<"(13 ("<<to_hex(fd.zone_index)<<" "<<to_hex(fd.istart+1)<<" "<<to_hex(fd.iend)<<" ";
		fs<<to_hex(fd.zone_type)<<" 2)(\n";

		for (int i=fd.istart; i<fd.iend; ++i){
			auto& ed = *g.vedges[iface++];
			int ip1 = ed.first()->id;
			int ip2 = ed.last()->id;
			int ic1 = ed.has_left_cell() ? ed.left.lock()->id : -1;
			int ic2 = ed.has_right_cell() ? ed.right.lock()->id : -1;
			fs<<to_hex(ip1+1)<<" "<<to_hex(ip2+1)<<" ";
			fs<<to_hex(ic1+1)<<" "<<to_hex(ic2+1)<<"\n";
		}
		fs<<"))\n";
	}
	//Periodic features
	int it = 0;
	for (int k=0; k<pd.size(); ++k){
		auto z1 = FaceData::zone_by_bindex(face_data, pd.b1[k]);
		auto z2 = FaceData::zone_by_bindex(face_data, pd.b2[k]);
		int sz = z1->n;
		fs<<"(18 ("<<to_hex(it+1)<<" "<<to_hex(it+sz)<<" ";
		fs<<to_hex(z1->zone_index)<<" "<<to_hex(z2->zone_index)<<")(\n";
		for (int i=0; i<sz; ++i){
			fs<<to_hex(periodic_relations[2*(it+i)] + 1)<<" "<<to_hex(periodic_relations[2*(it+i)+1] + 1)<<"\n";
		}
		fs<<"))\n";
		it+=sz;
	}

	//Boundary features
	fs<<"(45 (2 fluid fluid 1)())\n";
	for (auto& fd: face_data){
		fs<<"(45 ("<<fd.zone_index<<" "<<fd.zone_name<<" "<<fd.boundary_name<<" 1)())\n";
	}
	
	fs.close();
}

void hme::GridMSH(const GridData& g, std::string fn){
	return hme::GridMSH(g, fn, default_bfun, PeriodicData());
}

void hme::GridMSH(const GridData& g, std::string fn, PeriodicData pd){
	return hme::GridMSH(g, fn, default_bfun, pd);
}

void hme::GridMSH(const GridData& g, std::string fn, hme::BNamesFun bnames){
	return hme::GridMSH(g, fn, bnames, PeriodicData());
}
