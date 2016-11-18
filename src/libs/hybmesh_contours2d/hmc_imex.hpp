#ifndef HYBMESH_HMC_EXPORT_HPP
#define HYBMESH_HMC_EXPORT_HPP

#include "hmxmlreader.hpp"
#include "collections.hpp"
#include "containers.hpp"

namespace HMCont2D{ namespace Export{

struct EColWriter{
	//tp = "ascii", "bin", "floatbin"
	EColWriter(const ECollection& g,
			HMXML::ReaderA* writer,
			HMXML::Reader* subnode,
			std::string contname,
			std::string tp);

	void AddVertexData(std::string fieldname, const vector<double>& data, bool binary);
	void AddVertexData(std::string fieldname, const vector<float>& data, bool binary);
	void AddVertexData(std::string fieldname, const vector<char>& data, bool binary);
	void AddVertexData(std::string fieldname, const vector<int>& data, bool binary);

	void AddEdgeData(std::string fieldname, const vector<double>& data, bool binary);
	void AddEdgeData(std::string fieldname, const vector<float>& data, bool binary);
	void AddEdgeData(std::string fieldname, const vector<char>& data, bool binary);
	void AddEdgeData(std::string fieldname, const vector<int>& data, bool binary);

	void AddVertexData(std::string fieldname, const vector<vector<double>>& data, bool binary);
	void AddVertexData(std::string fieldname, const vector<vector<float>>& data, bool binary);
	void AddVertexData(std::string fieldname, const vector<vector<char>>& data, bool binary);
	void AddVertexData(std::string fieldname, const vector<vector<int>>& data, bool binary);

	void AddEdgeData(std::string fieldname, const vector<vector<double>>& data, bool binary);
	void AddEdgeData(std::string fieldname, const vector<vector<float>>& data, bool binary);
	void AddEdgeData(std::string fieldname, const vector<vector<char>>& data, bool binary);
	void AddEdgeData(std::string fieldname, const vector<vector<int>>& data, bool binary);

	template<class A> bool is_binary(){return false;}
private:
	virtual void data_changed(){}

	HMXML::ReaderA* pwriter;
	const ECollection* cont;
	HMXML::Reader cwriter, vwriter, ewriter;
	std::string __tp;
};
template<> bool EColWriter::is_binary<int>(){return __tp == "bin";}
template<> bool EColWriter::is_binary<char>(){return __tp == "bin";}
template<> bool EColWriter::is_binary<double>(){return __tp != "ascii";}
template<> bool EColWriter::is_binary<float>(){return __tp != "ascii";}

}}

namespace HMCont2D{namespace Import{
struct EColReader{
	EColReader(HMXML::ReaderA* preader, HMXML::Reader* subnode);
	std::unique_ptr<HMCont2D::Container<HMCont2D::ECollection>> result;

	struct TFieldInfo{
		TFieldInfo(HMXML::Reader& field);
		std::string name, type;
		int dim;
	};
	std::vector<TFieldInfo> edges_fields();
	std::vector<TFieldInfo> vertices_fields();

	template<class A> vector<A> read_vertices_field(std::string name);
	template<class A> vector<A> read_edges_field(std::string name);

	template<class A> vector<vector<A>> read_vertices_vecfield(std::string name);
	template<class A> vector<vector<A>> read_edges_vecfield(std::string name);

	int Ne, Nv;
private:
	void fill_result();
	HMXML::Reader* pgreader;
	HMXML::ReaderA* preader;

	template<class A>
	vector<A> read_field(HMXML::Reader& rd, int num);
	template<class A>
	vector<vector<A>> read_vecfield(HMXML::Reader& rd, int num);
};

template<class A>
vector<A> EColReader::read_vertices_field(std::string name){
	auto curnode = pgreader->find_by_path("VERTICES/FIELD[@name='"+name+"']", true);
	return read_field<A>(curnode, Nv);
}
template<class A>
vector<A> EColReader::read_edges_field(std::string name){
	auto curnode = pgreader->find_by_path("EDGES/FIELD[@name='"+name+"']", true);
	return read_field<A>(curnode, Ne);
}
template<class A>
vector<A> EColReader::read_field(HMXML::Reader& rd, int num){
	TFieldInfo info(rd);
	auto content = preader->read_num_content(rd, num);
	if (info.type == "int") return content.convert_data<int, A>();
	if (info.type == "float") return content.convert_data<float, A>();
	if (info.type == "double") return content.convert_data<double, A>();
	if (info.type == "char") return content.convert_data<char, A>();
	return vector<A>();
}

template<class A>
vector<vector<A>> EColReader::read_vertices_vecfield(std::string name){
	auto curnode = pgreader->find_by_path("VERTICES/FIELD[@name='"+name+"']", true);
	return read_vecfield<A>(curnode, Nv);
}
template<class A>
vector<vector<A>> EColReader::read_edges_vecfield(std::string name){
	auto curnode = pgreader->find_by_path("EDGES/FIELD[@name='"+name+"']", true);
	return read_vecfield<A>(curnode, Ne);
}
template<class A>
vector<vector<A>> EColReader::read_vecfield(HMXML::Reader& rd, int num){
	TFieldInfo info(rd);
	auto content = preader->read_num_content(rd, num);
	if (info.type == "int") return content.convert_vdata<int, A>();
	if (info.type == "float") return content.convert_vdata<float, A>();
	if (info.type == "double") return content.convert_vdata<double, A>();
	if (info.type == "char") return content.convert_vdata<char, A>();
	return vector<vector<A>>();
}

}}


#endif
