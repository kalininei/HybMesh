#include "hmc_imex.hpp"

using namespace HMCont2D;
using namespace HMXML;

namespace{
template<class A>
void add_field(Reader& subnode, std::string fieldname, const vector<A>& data, bool binary, ReaderA& writer){
	Reader wr = subnode.find_by_path("FIELD[@name='"+fieldname+"']");
	if (wr) wr.unlink_node();
	wr = subnode.new_child("FIELD");
	wr.new_attribute("name", fieldname);
	writer.set_num_content(data, wr, binary);
}
}

Export::EColWriter::EColWriter(const ECollection& c,
		HMXML::ReaderA* writer,
		HMXML::Reader* subnode,
		std::string contname,
		std::string tp): pwriter(writer), cont(&c){
	__tp = tp;
	
	//supplementary data
	auto ap = c.all_points();
	vector<double> pcoords(ap.size()*2);
	for (int i=0; i<ap.size(); ++i){
		pcoords[2*i] = ap[i]->x;
		pcoords[2*i+1] = ap[i]->y;
	}
	vector<int> edgeconnect; edgeconnect.reserve(c.size()*2);
	auto _indexer = aa::ptr_container_indexer(ap);
	_indexer.convert();
	for (auto e: *cont){
		edgeconnect.push_back(_indexer.index(e->pstart));
		edgeconnect.push_back(_indexer.index(e->pend));
	}
	_indexer.restore();

	//create xml structure
	cwriter = subnode->new_child("CONTOUR2D");
	cwriter.new_attribute("name", contname);

	//inforamtion
	cwriter.new_child("N_VERTICES").set_content(std::to_string(pcoords.size()/2));
	cwriter.new_child("N_EDGES").set_content(std::to_string(cont->size()));

	//vertices
	vwriter = cwriter.new_child("VERTICES");
	auto coordswriter = vwriter.new_child("COORDS");
	writer->set_num_content(pcoords, coordswriter, is_binary<double>());

	//edges
	ewriter = cwriter.new_child("EDGES");
	auto vconnectwriter = ewriter.new_child("VERT_CONNECT");
	writer->set_num_content(edgeconnect, vconnectwriter, is_binary<int>());
}


void Export::EColWriter::AddVertexData(std::string fieldname, const vector<double>& data, bool binary){
	data_changed();
	add_field<double>(vwriter, fieldname, data, binary, *pwriter);
}
void Export::EColWriter::AddVertexData(std::string fieldname, const vector<float>& data, bool binary){
	data_changed();
	add_field<float>(vwriter, fieldname, data, binary, *pwriter);
}
void Export::EColWriter::AddVertexData(std::string fieldname, const vector<char>& data, bool binary){
	data_changed();
	add_field<char>(vwriter, fieldname, data, binary, *pwriter);
}
void Export::EColWriter::AddVertexData(std::string fieldname, const vector<int>& data, bool binary){
	data_changed();
	add_field<int>(vwriter, fieldname, data, binary, *pwriter);
}

void Export::EColWriter::AddEdgeData(std::string fieldname, const vector<double>& data, bool binary){
	data_changed();
	add_field<double>(ewriter, fieldname, data, binary, *pwriter);
}
void Export::EColWriter::AddEdgeData(std::string fieldname, const vector<float>& data, bool binary){
	data_changed();
	add_field<float>(ewriter, fieldname, data, binary, *pwriter);
}
void Export::EColWriter::AddEdgeData(std::string fieldname, const vector<char>& data, bool binary){
	data_changed();
	add_field<char>(ewriter, fieldname, data, binary, *pwriter);
}
void Export::EColWriter::AddEdgeData(std::string fieldname, const vector<int>& data, bool binary){
	data_changed();
	add_field<int>(ewriter, fieldname, data, binary, *pwriter);
}
void Export::EColWriter::AddVertexData(std::string fieldname, const vector<vector<double>>& data, bool binary){
	data_changed();
	add_field<vector<double>>(vwriter, fieldname, data, binary, *pwriter);
}
void Export::EColWriter::AddVertexData(std::string fieldname, const vector<vector<float>>& data, bool binary){
	data_changed();
	add_field<vector<float>>(vwriter, fieldname, data, binary, *pwriter);
}
void Export::EColWriter::AddVertexData(std::string fieldname, const vector<vector<char>>& data, bool binary){
	data_changed();
	add_field<vector<char>>(vwriter, fieldname, data, binary, *pwriter);
}
void Export::EColWriter::AddVertexData(std::string fieldname, const vector<vector<int>>& data, bool binary){
	data_changed();
	add_field<vector<int>>(vwriter, fieldname, data, binary, *pwriter);
}

void Export::EColWriter::AddEdgeData(std::string fieldname, const vector<vector<double>>& data, bool binary){
	data_changed();
	add_field<vector<double>>(ewriter, fieldname, data, binary, *pwriter);
}
void Export::EColWriter::AddEdgeData(std::string fieldname, const vector<vector<float>>& data, bool binary){
	data_changed();
	add_field<vector<float>>(ewriter, fieldname, data, binary, *pwriter);
}
void Export::EColWriter::AddEdgeData(std::string fieldname, const vector<vector<char>>& data, bool binary){
	data_changed();
	add_field<vector<char>>(ewriter, fieldname, data, binary, *pwriter);
}
void Export::EColWriter::AddEdgeData(std::string fieldname, const vector<vector<int>>& data, bool binary){
	data_changed();
	add_field<vector<int>>(ewriter, fieldname, data, binary, *pwriter);
}


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
	
	result.reset(new Container<ECollection>());
	//add points
	for (int i=0; i<Nv; ++i){
		result->pdata.add_value(Point(vert[2*i], vert[2*i+1]));
	}
	//add edges
	for (int i=0; i<Ne; ++i){
		auto p1 = result->pdata.pvalue(edgevert[2*i]);
		auto p2 = result->pdata.pvalue(edgevert[2*i+1]);
		result->add_value(Edge(p1, p2));
	}
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
