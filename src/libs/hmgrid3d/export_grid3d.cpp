#include "export_grid3d.hpp"
#include <unordered_map>
#include <fstream>
#include <sstream>
#include "debug_grid3d.hpp"
using namespace HMGrid3D;
namespace hme = HMGrid3D::Export;

std::tuple<
	std::vector<double>,   //points coordinates in plain array
	std::vector<int>       //????
>
hme::Serialize(const HMGrid3D::Grid& g){
	_THROW_NOT_IMP_;
}

namespace{

std::vector<int>
face_to_point_representation(const HMGrid3D::Face& f, std::map<HMGrid3D::Vertex*, int>& pind){
	std::vector<int> ret;
	for (auto p: f.sorted_vertices()) ret.push_back(pind[p.get()]);
	return ret;
}

std::vector<std::vector<int>>
cell_to_point_representation(const HMGrid3D::Cell& c, std::map<HMGrid3D::Vertex*, int>& pind){
	std::vector<std::vector<int>> ret;
	for (auto f: c.faces){
		ret.push_back(face_to_point_representation(*f, pind));
		//all faces have the cell at their left side
		if (f->left.get() != &c) std::reverse(ret.back().begin(), ret.back().end());
	}
	return ret;
}

struct vtkcell_expression{
	vector<int> pts;
	int celltype;

	int wsize() const { return 1+pts.size(); }
	std::string to_string() const {
		std::ostringstream out;
		out<<pts.size();
		for (auto i:pts) out<<" "<<i;
		return out.str();
	}
	std::string stype() const { return std::to_string(celltype); }

	bool _false_return(){
		pts.clear();
		celltype = 0;
		return false;
	}

	static int find_opposite(std::vector<std::vector<int>>& data, int cind, int pind){
		//find point opposite to data[cind][pind]
		int p_oppose = -1;
		for (int i=0; i<data.size(); ++i) if (i!=cind){
			auto fnd = std::find(data[i].begin(), data[i].end(), data[cind][pind]);
			if (fnd == data[i].end()) continue;
			auto fndnext = std::next(fnd);
			if (fndnext == data[i].end()) fndnext = data[i].begin();
			if (std::find(data[cind].begin(), data[cind].end(), *fndnext) ==
					data[cind].end()){
				p_oppose = *fndnext;
				break;
			}
			auto fndprev = (fnd == data[i].begin()) ? --data[i].end() : std::prev(fnd);
			if (std::find(data[cind].begin(), data[cind].end(), *fndprev) ==
					data[cind].end()){
				p_oppose = *fndprev;
				break;
			}
		}
		return p_oppose;
	}
	static int is_opposite(std::vector<std::vector<int>>& data, int i1, int i2){
		//returns -1 or index of data[i2] which matches data[i1][0]
		int p_oppose = find_opposite(data, i1, 0);
		if (p_oppose < 0) return -1;
		auto fnd = std::find(data[i2].begin(), data[i2].end(), p_oppose);
		if (fnd != data[i2].end()) return fnd - data[i2].begin();
		else return -1;
	}
	static std::pair<int, int> get_opposite(std::vector<std::vector<int>>& data, int ic){
		// returns index of opposite face, number of edge within data[first]
		// array which mathces data[ic][0] vertex. (-1, -1) if no oppsite was found
		std::pair<int, int> bad_ret(-1, -1);
		//find vertex opposite to data[ic][0] as one that is connected to data[ic][0]
		//but not lieing in data[ic]
		int p_oppose = find_opposite(data, ic, 0);
		if (p_oppose<0) return bad_ret;
		//find opposite face as face not containing any data[ic] nodes
		//but containing p_oppose
		for (auto i=0; i<data.size(); ++i) if (i != ic){
			//same number of points
			if (data[i].size() != data[ic].size()) continue;
			//no dublicate points
			std::unordered_set<int> dub(data[i].begin(), data[i].end());
			dub.insert(data[ic].begin(), data[ic].end());
			if (dub.size() != data[i].size() + data[ic].size()) continue;
			//find p_oppose
			auto fnd = std::find(data[i].begin(), data[i].end(), p_oppose);
			if (fnd == data[i].end()) continue;
			return std::make_pair(i, fnd-data[i].begin());
		}
		return bad_ret;
	}
	bool try_tetrahedron(std::vector<std::vector<int>>& data){
		celltype = 10;
		if (data.size() != 4) return _false_return();
		for (auto& d: data) if (d.size() != 3) return _false_return();
		for (auto& d: data) pts.insert(pts.end(), d.begin(), d.end()); 
		pts = aa::no_dublicates(pts);
		if (pts.size() == 4) return true;
		else return _false_return();
	}
	bool try_hexahedron(std::vector<std::vector<int>>& data){
		celltype = 12;
		if (data.size() != 6) return _false_return();
		for (auto& d: data) if (d.size() != 4) return _false_return();
		//find opposite face
		std::pair<int, int> opcell = get_opposite(data, 0);
		if (opcell.first < 0) return _false_return();
		std::vector<int> dlower = data[0];
		std::vector<int> dupper = data[opcell.first];
		std::rotate(dupper.begin(), dupper.begin()+opcell.second, dupper.end());
		pts = {dlower[0], dlower[3], dlower[2], dlower[1],
			dupper[0], dupper[1], dupper[2], dupper[3]};
		return true;
	}
	bool try_wedge(std::vector<std::vector<int>>& data){
		celltype=13; 
		if (data.size() != 5) return _false_return();
		int f1=-1, f2=-1;
		for (int i=0; i<data.size(); ++i){
			auto& c = data[i];
			if (c.size() <3 || c.size() > 4) return _false_return();
			if (c.size() == 4) continue;
			if (f1<0) f1 = i;
			else if (f2<0) f2 = i;
			else return _false_return();
		}
		if (f2<0) return _false_return();
		auto oi = is_opposite(data, f1, f2);
		if (oi < 0) return _false_return();
		auto dlower = data[f1], dupper = data[f2];
		std::rotate(dupper.begin(), dupper.begin() + oi, dupper.end());
		pts = {dlower[0], dlower[1], dlower[2], dupper[0], dupper[2], dupper[1]};
		return true;
	}
};

vtkcell_expression get_vtkcell(const HMGrid3D::Cell& c, std::map<HMGrid3D::Vertex*, int>& pind){
	vtkcell_expression ret;
	auto cint = cell_to_point_representation(c, pind);
	if (ret.try_tetrahedron(cint)) return ret;
	if (ret.try_hexahedron(cint)) return ret;
	if (ret.try_hexahedron(cint)) return ret;
	if (ret.try_wedge(cint)) return ret;
	std::string s("Cannot write 3D cell with ");
	s += std::to_string(c.n_faces());
	s += " faces to vtk format";
	throw std::runtime_error(s.c_str());
}

}

void hme::GridVTK(const HMGrid3D::Grid& g, const char* fn){
	//grid data by indices
	auto pts = g.allvertices();
	auto cells = g.allcells();
	std::map<HMGrid3D::Vertex*, int> pindall;
	int i=0; for (auto p: pts) pindall[p.get()] = i++;

	//write to file:
	//header
	std::ofstream fs(fn);
	fs<<"# vtk DataFile Version 3.0"<<std::endl;
	fs<<"3D Grid"<<std::endl;
	fs<<"ASCII"<<std::endl;

	//Points
	fs<<"DATASET UNSTRUCTURED_GRID"<<std::endl;
	fs<<"POINTS "<<pts.size()<< " float"<<std::endl;
	for (auto p: pts) fs<<p->x<<" "<<p->y<<" "<<p->z<<std::endl;

	//Cells
	vector< vtkcell_expression > vtkcell;
	for (auto c: cells) vtkcell.push_back(get_vtkcell(*c, pindall));
	int nffull = std::accumulate(vtkcell.begin(), vtkcell.end(), 0, 
			[](int s, const vtkcell_expression& v){ return s + v.wsize(); });
	fs<<"CELLS  "<<vtkcell.size()<<"   "<<nffull<<std::endl;
	for (auto& f: vtkcell) fs<<f.to_string()<<std::endl;
	fs<<"CELL_TYPES  "<<vtkcell.size()<<std::endl;
	for (auto& f: vtkcell) fs<<f.stype()<<std::endl;

	fs.close();
}

void hme::BoundaryVTK(const HMGrid3D::Grid& g, const char* fn){
	//grid data by indices
	std::map<HMGrid3D::Vertex*, int> pindall;
	std::map<HMGrid3D::Face*, int> findall;
	std::map<HMGrid3D::Cell*, int> cindall;
	int i=0; for (auto p: g.allvertices()) pindall[p.get()] = i++;
	    i=0; for (auto f: g.allfaces()) if (f->boundary_type >= 0) findall[f.get()] = i++;
	    i=0; for (auto c: g.allcells()) cindall[c.get()] = i++;

	//boundary points and faces
	vector<Face*> bfaces;
	for (auto f: g.allfaces()) if (f->is_boundary()) bfaces.push_back(f.get());
	vector<Vertex*> bpoints;
	std::set<Vertex*> usedp;
	for (auto f: bfaces){
		for (auto p: f->allvertices()){
			if (usedp.emplace(p.get()).second) bpoints.push_back(p.get());
		}
	}
	//point in bpoints -> its index in bpoints
	std::map<HMGrid3D::Vertex*, int> pindlocal;
	i=0; for (auto p: bpoints) pindlocal[p] = i++;
	//fill faces array by vertex indices
	vector<vector<int>> outfaces;
	for (auto f: bfaces) outfaces.push_back(face_to_point_representation(*f, pindlocal));

	//fill additional data arrays
	vector<int> f_btypes, f_gindices, p_gindices, f_cellindices;
	for (auto p: bpoints) p_gindices.push_back(pindall[p]);
	for (auto f: bfaces){
		f_gindices.push_back(findall[f]);
		f_btypes.push_back(f->boundary_type);
		auto cell = (f->left) ? f->left : f->right;
		f_cellindices.push_back(cindall[cell.get()]);
	}

	//write to file:
	//header
	std::ofstream fs(fn);
	fs<<"# vtk DataFile Version 3.0"<<std::endl;
	fs<<"Boundary for 3D Grid"<<std::endl;
	fs<<"ASCII"<<std::endl;

	//Points
	fs<<"DATASET UNSTRUCTURED_GRID"<<std::endl;
	fs<<"POINTS "<<bpoints.size()<< " float"<<std::endl;
	for (auto p: bpoints) fs<<p->x<<" "<<p->y<<" "<<p->z<<std::endl;

	//Faces
	int nffull = std::accumulate(outfaces.begin(), outfaces.end(), 0, 
			[](int s, const vector<int>& v){ return s + v.size(); });
	fs<<"CELLS  "<<outfaces.size()<<"   "<<nffull + outfaces.size()<<std::endl;
	for (auto& f: outfaces){
		fs<<f.size();
		for (auto ip: f) fs<<" "<<ip; fs<<std::endl;
	}
	fs<<"CELL_TYPES  "<<outfaces.size()<<std::endl;
	for (int i=0;i<outfaces.size();++i) fs<<7<<" "; fs<<std::endl;

	//additional info
	auto wintarray = [&](const vector<int>& dt, const char* nm){
		fs<<"SCALARS "<<nm<<" int 1"<<std::endl;
		fs<<"LOOKUP_TABLE default"<<std::endl;
		for (auto d: dt) fs<<d<<" "; fs<<std::endl;
	};

	//faces boundary types
	fs<<"CELL_DATA"<<" "<<f_btypes.size()<<std::endl;
	wintarray(f_btypes, "boundaries");
	//faces global indices
	wintarray(f_gindices, "face_global_indices");
	//face adjacent cell
	wintarray(f_cellindices, "adjacent_cell_indices");
	//point global indices
	fs<<"POINT_DATA"<<" "<<p_gindices.size()<<std::endl;
	wintarray(p_gindices, "vertex_global_indices");

	fs.close();
}

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
	int nf = c.n_faces();
	int ne = c.n_edges();
	int nv = c.n_vertices();
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

//save to fluent main function
void gridmsh(const HMGrid3D::Grid& g, const char* fn,
		std::function<std::string(int)> btype_name,
		std::map<Face*, Face*> pfaces=std::map<Face*, Face*>()){
	if (g.n_cells() == 0) throw std::runtime_error("Exporting blank grid");
	//Zones:
	//1    - verticies default
	//2    - fluid for cells
	//3    - default interior
	//4..N - bc's faces

	//Needed data
	auto vert = g.allvertices();
	auto cells = g.allcells();
	std::map<int, ShpVector<Face>> facezones;  //faces by zones id
	{
		//interior zone has always be there
		facezones.emplace(std::numeric_limits<int>::min(), ShpVector<Face>());
		for (auto f: g.allfaces()){
			int tp = (f->is_boundary()) ? f->boundary_type : std::numeric_limits<int>::min();
			auto emp = facezones.emplace(tp, ShpVector<Face>());
			emp.first->second.push_back(f);
		}
	}
	int nfaces=0;                 //total number of faces
	{
		for (auto& fz: facezones) nfaces+=fz.second.size();
	}
	std::vector<std::pair<int, int>> face_adjacents; //left/right cell indicies for each cell
	{
		std::unordered_map<HMGrid3D::Cell*, int> cell_ind;
		int i=0;
		for (auto c: cells) cell_ind.emplace(c.get(), i++);
		for (auto& fz: facezones) for (auto f: fz.second){
			face_adjacents.push_back(std::make_pair(-1, -1));
			if (f->right) face_adjacents.back().first = cell_ind[f->right.get()];
			if (f->left) face_adjacents.back().second = cell_ind[f->left.get()];
		}
	}
	std::vector<std::vector<int>> face_vertices;
	{
		std::unordered_map<Vertex*, int> vert_ind;
		int i=0;
		for (auto p: vert) vert_ind.emplace(p.get(), i++); 
		for (auto& fz: facezones) for (auto f: fz.second){
			face_vertices.push_back(vector<int>());
			for (auto p: f->sorted_vertices()){
				face_vertices.back().push_back(vert_ind[p.get()]);
			}
		}
	}
	std::vector<char> ctypes;   //cell types array
	{
		for (auto c: cells) ctypes.push_back(msh_cell_type(*c));
	}
	char cell_common_type;     //common cell types or '0' if they differ
	{
		cell_common_type = ctypes[0];
		for (auto c: ctypes){
			if (c != cell_common_type){
				cell_common_type = '0';
				break;
			}
		}
	}
	std::vector<char> ftypes;   //face types for each zone
	{
		for (auto& fz: facezones){
			for (auto f: fz.second) ftypes.push_back(msh_face_type(*f));
		}
	}
	std::vector<char> facezones_common_type;  //common face type or '0' for each zone
	{
		int iface = 0;
		for (auto& fz: facezones){
			facezones_common_type.push_back(ftypes[iface]);
			for (int i=0; i<fz.second.size(); ++i){
				if (ftypes[iface+i] != facezones_common_type.back()){
					facezones_common_type.back() = '0';
					break;
				}
			}
			iface += fz.second.size();
		}
	}
	std::map<
		std::pair<int, int>,
		std::vector<std::pair<int, int>>
	> periodic_data;     //<periodic zonetype, shadow zonetype> -> vector of periodic face, shadow face>
	{
		auto perface_index = [&](Face* f){
			int ret = 0;
			auto fnd1 = facezones.find(f->boundary_type);
			for (auto it=facezones.begin(); it!= fnd1; ++it){
				ret += it->second.size();
			}
			auto fnd3 = aa::shp_find(fnd1->second.begin(), fnd1->second.end(), f);
			return ret + fnd3 - (fnd1->second.begin());
		};
		std::map<int, int> btype_zone;
		int zt = 4;
		for (auto it=++facezones.begin(); it!=facezones.end(); ++it){
			btype_zone[it->first] = zt++;
		}
		for (auto& ff: pfaces){
			int perzone = btype_zone[ff.first->boundary_type];
			int shazone = btype_zone[ff.second->boundary_type];
			int perindex = perface_index(ff.first);
			int shaindex = perface_index(ff.second);
			auto emp = periodic_data.emplace(
					std::make_pair(perzone, shazone),
					std::vector<std::pair<int, int>>());
			emp.first->second.emplace_back(perindex, shaindex);
		}
	}
	std::map<int, std::string> facezones_btype; //2-interior, 3-wall, 8-periodic-shadow, 12-periodic
	std::map<int, std::string> facezones_bname;
	{
		std::set<int> periodic_zones, pershadow_zones;
		for (auto& d: periodic_data){
			periodic_zones.insert(d.first.first); 
			pershadow_zones.insert(d.first.second); 
		}
		facezones_btype[3] = to_hex(2);
		facezones_bname[3] = "interior";
		int zt = 4;
		for (auto it = ++facezones.begin(); it!=facezones.end(); ++it){
			if (periodic_zones.find(it->first) != periodic_zones.end()){
				facezones_btype[zt] = to_hex(12);
				facezones_bname[zt] = "periodic";
			} else if (pershadow_zones.find(it->first) != pershadow_zones.end()){
				facezones_btype[zt] = to_hex(8);
				facezones_bname[zt] = "periodic-shadow";
			} else {
				facezones_btype[zt] = to_hex(3);
				facezones_bname[zt] = "wall";
			}
			++zt;
		}
	}
	
	//=========== Write to file
	std::ofstream fs(fn);
	fs.precision(17);
	//header
	fs<<"(0 \"HybMesh to Fluent File\")\n(0 \"Dimensions\")\n(2 3)\n";

	//Vertices: Zone 1
	fs<<"(10 (0 1 "<<to_hex(vert.size())<<" 0 3))\n";
	fs<<"(10 (1 1 "<<to_hex(vert.size())<<" 1 3)(\n";
	for (auto p: vert){
		fs<<p->x<<" "<<p->y<<" "<<p->z<<"\n";
	}
	fs<<"))\n";

	//Cells: Zone 2
	fs<<"(12 (0 1 "<<to_hex(cells.size())<<" 0))\n";
	fs<<"(12 (2 1 "<<to_hex(cells.size())<<" 1 "<<cell_common_type<<")";
	if (cell_common_type != '0') fs<<")\n";
	else {
		fs<<"(\n";
		for (auto s: ctypes) fs<<s<<" ";
		fs<<"\n))\n";
	}

	//Faces: Zones 3+it
	fs<<"(13 (0 1 "<<to_hex(nfaces)<<" 0))\n";
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
			fs<<"(18 ("<<to_hex(it+1)<<" "<<to_hex(it+sz+1)<<" ";
			fs<<to_hex(pd.first.first)<<" "<<to_hex(pd.first.second)<<")(\n";
			for (auto& v: pd.second){
				fs<<to_hex(v.first+1)<<" "<<to_hex(v.second+1)<<"\n";
			}
			fs<<")\n";
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

HMGrid3D::Grid hme::PeriodicData::assemble(const HMGrid3D::Grid& g,
		std::map<Face*, Face*>& outmap){
	if (size() == 0) return g;
	_THROW_NOT_IMP_;
}

void hme::GridMSH(const HMGrid3D::Grid& g, const char* fn,
		std::function<std::string(int)> btype_name,
		PeriodicData periodic){
	if (periodic.size() == 0) return gridmsh(g, fn, btype_name);
	else{
		std::map<Face*, Face*> periodic_cells;
		HMGrid3D::Grid gp = periodic.assemble(g, periodic_cells);
		return gridmsh(g, fn, btype_name, periodic_cells);
	}
}

void hme::GridMSH(const HMGrid3D::Grid& g, const char* fn){
	return gridmsh(g, fn, [](int i)->std::string{
			return std::string("boundary") + std::to_string(i);
		;});
}

void hme::GridMSH(const HMGrid3D::Grid& g, const char* fn,
		std::function<std::string(int)> btype_name){
	return gridmsh(g, fn, btype_name);
}
