#include "import2d_hm.hpp"
#include "constructor.hpp"

using namespace HM2D;
using namespace HMXML;

Import::EColReader::EColReader(HMXML::ReaderA* reader, HMXML::Reader* subnode):preader(reader), pgreader(subnode){
	fill_result();
}
void Import::EColReader::fill_result(){
	//read dimensions
	pgreader->value_int("N_VERTICES", Nv, true);
	pgreader->value_int("N_EDGES", Ne, true);
	//read vertices
	Reader tmp = pgreader->find_by_path("VERTICES/COORDS", true);
	vector<double> vert = preader->read_num_content(tmp, 2*Nv).vec<double>();
	//read edges
	tmp = pgreader->find_by_path("EDGES/VERT_CONNECT", true);
	vector<int> edgevert = preader->read_num_content(tmp, 2*Ne).vec<int>();
	for (auto i: edgevert) if (i>=Nv || i<0) throw std::runtime_error("Edge-Vertex connectivity contains illegal vertex index");
	
	result.reset(new EdgeData);
	*result = HM2D::ECol::Constructor::FromRaw(Nv, Ne, &vert[0], &edgevert[0]);
}
std::vector<Import::EColReader::TFieldInfo> Import::EColReader::edges_fields(){
	std::vector<EColReader::TFieldInfo> ret;
	auto fnd = pgreader->findall_by_path("EDGES/FIELD");
	for (auto& f: fnd) ret.push_back(TFieldInfo(f));
	return ret;
}
std::vector<Import::EColReader::TFieldInfo> Import::EColReader::vertices_fields(){
	std::vector<EColReader::TFieldInfo> ret;
	auto fnd = pgreader->findall_by_path("VERTICES/FIELD");
	for (auto& f: fnd) ret.push_back(TFieldInfo(f));
	return ret;
}

Import::EColReader::TFieldInfo::TFieldInfo(HMXML::Reader& field){
	name = field.attribute(".", "name");
	type = field.attribute(".", "type");
	std::string sdim = field.attribute(".", "dim");

	if (name=="" || type=="") throw std::runtime_error(
			"not enough attributes in FIELD node");
	dim = 1;
	if (sdim=="variable") dim=-1;
	if (sdim!="") dim=atoi(sdim.c_str());
}

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

GridData Import::GridFromTabs(const vector<double> vert, const vector<int>& edgevert,
		const vector<int>& edgecell){
	int Nv = vert.size()/2;
	int Nc = *std::max_element(edgecell.begin(), edgecell.end())+1;
	int Ne = edgevert.size()/2;
	GridData ng;
	for (int i=0; i<Nv; ++i){
		ng.vvert.push_back(std::make_shared<Vertex>(vert[2*i], vert[2*i+1]));
	}
	for (int i=0; i<Ne; ++i){
		int p1 = edgevert[2*i], p2 = edgevert[2*i+1];
		ng.vedges.push_back(std::make_shared<Edge>(ng.vvert[p1], ng.vvert[p2]));
	}
	for (int i=0; i<Nc; ++i){ ng.vcells.push_back(std::make_shared<Cell>()); }
	for (int i=0; i<Ne; ++i){
		int c1 = edgecell[2*i], c2 = edgecell[2*i+1];
		if (c1>=0){
			ng.vcells[c1]->edges.push_back(ng.vedges[i]);
			ng.vedges[i]->left = ng.vcells[c1];
		}
		if (c2>=0){
			ng.vcells[c2]->edges.push_back(ng.vedges[i]);
			ng.vedges[i]->right = ng.vcells[c2];
		}
	}
	//sort cells->edges tables
	aa::enumerate_ids_pvec(ng.vvert);
	for (int i=0; i<Nc; ++i){
		auto& c = ng.vcells[i];
		vector<int> p1(c->edges.size());
		vector<int> p2(c->edges.size());
		for (int j=0; j<c->edges.size(); ++j){
			p1[j] = c->edges[j]->pfirst()->id;
			p2[j] = c->edges[j]->plast()->id;
			if (c->edges[j]->left.lock()!=c){
				std::swap(p1[j], p2[j]);
			}
		}
		for (int j=1; j<c->edges.size(); ++j){
			int p1s = p2[j-1];
			int jfnd = std::find(p1.begin()+j, p1.end(), p1s) - p1.begin();
			if (jfnd != j){
				if (jfnd >= p1.size()) throw std::runtime_error("Invalid grid geometry");
				std::swap(c->edges[j], c->edges[jfnd]);
				std::swap(p1[j], p1[jfnd]);
				std::swap(p2[j], p2[jfnd]);
			}
		}
	}
	return ng;
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

	//constructing a grid
	GridData ng = Import::GridFromTabs(vert, edgevert, edgecell);
	result.reset(new GridData(std::move(ng)));
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
