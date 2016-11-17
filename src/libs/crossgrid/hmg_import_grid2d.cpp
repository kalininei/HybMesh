#include "hmg_import_grid2d.hpp"
using namespace GGeom;
using namespace HMXML;

Import::WOwner::WOwner(std::string fn, std::string name){
	reader.reset(new ReaderA(fn, "</HybMeshData>"));
	//find a node
	if (name==""){
		greader = reader->find_by_path("GRID2D");
	} else {
		greader = reader->find_by_path("GRID2D[@name='"+name+"']");
	}
	if (!greader) throw std::runtime_error("grid "+name+" was not found");
}

Import::GridReader::GridReader(HMXML::ReaderA* reader, HMXML::Reader* subnode):preader(reader), pgreader(subnode){
	fill_result();
}
void Import::GridReader::fill_result(){
	//read dimensions
	pgreader->value_int("N_VERTICES", Nv, true);
	pgreader->value_int("N_EDGES", Ne, true);
	pgreader->value_int("N_CELLS", Nc, true);
	//read vertices
	Reader tmp = pgreader->find_by_path("VERTICES/COORDS", true);
	vector<double> vert = preader->read_num_content(tmp, 2*Nv).vec<double>();
	//read edges
	tmp = pgreader->find_by_path("EDGES/VERT_CONNECT", true);
	vector<int> edgevert = preader->read_num_content(tmp, 2*Ne).vec<int>();
	tmp = pgreader->find_by_path("EDGES/CELL_CONNECT", true);
	vector<int> edgecell = preader->read_num_content(tmp, 2*Ne).vec<int>();
	for (auto i: edgevert) if (i>=Nv || i<0) throw std::runtime_error("Edge-Vertex connectivity contains illegal vertex index");
	for (auto i: edgecell) if (i>=Nc) throw std::runtime_error("Edge-Cell connectivity contains illegal cell index");
	//construct cells_edges table
	std::vector<std::vector<int>> cells_edges(Nc);
	for (int i=0; i<edgecell.size()/2; ++i){
		if (edgecell[2*i]>=0) cells_edges[edgecell[2*i]].push_back(i);
		if (edgecell[2*i+1]>=0) cells_edges[edgecell[2*i+1]].push_back(i);
	}
	//renumber cells_edges table
	for (int i=0; i<cells_edges.size(); ++i){
		if (cells_edges.size()==0) continue;
		for (int j=1; j<cells_edges[i].size(); ++j){
			int e0 = cells_edges[i][j-1];
			int pcur = (edgecell[2*e0]==i)?edgevert[2*e0+1]:edgevert[2*e0];
			int ind, e1;
			for (ind=j; ind<cells_edges[i].size(); ++ind){
				e1 = cells_edges[i][ind];
				int pfirst = (edgecell[2*e1]==i)?edgevert[2*e1]:edgevert[2*e1+1];
				if (pfirst == pcur) break;
			}
			std::swap(cells_edges[i][j], cells_edges[i][ind]);
		}
	}
	
	//construct cell->vertices connectivity
	int totsize = 0;
	for (int i=0; i<cells_edges.size(); ++i) totsize+=cells_edges[i].size();
	std::vector<int> cells_vertices(totsize+Nc);
	auto it=cells_vertices.begin();
	for (int i=0; i<Nc; ++i){
		*it++ = cells_edges[i].size();
		for (int j=0; j<cells_edges[i].size(); ++j){
			int e0 = cells_edges[i][j];
			*it++ = (edgecell[2*e0]==i)?edgevert[2*e0]:edgevert[2*e0+1];
		}
	}
	result.reset(new GridGeom(Nv, Nc, &vert[0], &cells_vertices[0]));
}
std::vector<Import::GridReader::TFieldInfo> Import::GridReader::edges_fields(){
	std::vector<GridReader::TFieldInfo> ret;
	auto fnd = pgreader->findall_by_path("EDGES/FIELD");
	for (auto& f: fnd) ret.push_back(TFieldInfo(f));
	return ret;
}
std::vector<Import::GridReader::TFieldInfo> Import::GridReader::vertices_fields(){
	std::vector<GridReader::TFieldInfo> ret;
	auto fnd = pgreader->findall_by_path("VERTICES/FIELD");
	for (auto& f: fnd) ret.push_back(TFieldInfo(f));
	return ret;
}
std::vector<Import::GridReader::TFieldInfo> Import::GridReader::cells_fields(){
	std::vector<GridReader::TFieldInfo> ret;
	auto fnd = pgreader->findall_by_path("CELLS/FIELD");
	for (auto& f: fnd) ret.push_back(TFieldInfo(f));
	return ret;
}

Import::GridReader::TFieldInfo::TFieldInfo(HMXML::Reader& field){
	name = field.attribute(".", "name");
	type = field.attribute(".", "type");
	std::string sdim = field.attribute(".", "dim");

	if (name=="" || type=="") throw std::runtime_error(
			"not enough attributes in FIELD node");
	dim = 1;
	if (sdim=="variable") dim=-1;
	if (sdim!="") dim=atoi(sdim.c_str());
}
