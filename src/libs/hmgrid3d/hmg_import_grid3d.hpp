#ifndef HMG_IMPORT_GRID3D_HPP
#define HMG_IMPORT_GRID3D_HPP
#include "serialize_grid3d.hpp"
#include "hmxmlreader.hpp"
#include "hmcallback.hpp"

namespace HMGrid3D{ namespace Import{
struct GridReader;
struct TReadHMG;

struct GridReader{
	GridReader(HMXML::ReaderA* reader, HMXML::Reader* subnode):
		preader(reader), pgreader(subnode){}
	//holds data
	std::unique_ptr<HMGrid3D::SGrid> result;

	struct TFieldInfo{
		TFieldInfo(HMXML::Reader& field);
		std::string name, type;
		int dim;
	};
	std::vector<TFieldInfo> edges_fields();
	std::vector<TFieldInfo> vertices_fields();
	std::vector<TFieldInfo> cells_fields();

	template<class A> vector<A> read_vertices_field(std::string name);
	template<class A> vector<A> read_edges_field(std::string name);
	template<class A> vector<A> read_faces_field(std::string name);
	template<class A> vector<A> read_cells_field(std::string name);

	template<class A> vector<vector<A>> read_vertices_vecfield(std::string name);
	template<class A> vector<vector<A>> read_edges_vecfield(std::string name);
	template<class A> vector<vector<A>> read_faces_vecfield(std::string name);
	template<class A> vector<vector<A>> read_cells_vecfield(std::string name);

	int Nc, Ne, Nv, Nf;
private:
	friend struct TReadHMG;

	void fill_result();
	HMXML::Reader* pgreader;
	HMXML::ReaderA* preader;

	template<class A>
	vector<A> read_field(HMXML::Reader& rd, int num);
	template<class A>
	vector<vector<A>> read_vecfield(HMXML::Reader& rd, int num);
};

template<class A>
vector<A> GridReader::read_vertices_field(std::string name){
	auto curnode = pgreader->find_by_path("VERTICES/FIELD[@name='"+name+"']", true);
	return read_field<A>(curnode, Nv);
}
template<class A>
vector<A> GridReader::read_edges_field(std::string name){
	auto curnode = pgreader->find_by_path("EDGES/FIELD[@name='"+name+"']", true);
	return read_field<A>(curnode, Ne);
}
template<class A>
vector<A> GridReader::read_faces_field(std::string name){
	auto curnode = pgreader->find_by_path("FACES/FIELD[@name='"+name+"']", true);
	return read_field<A>(curnode, Nc);
}
template<class A>
vector<A> GridReader::read_cells_field(std::string name){
	auto curnode = pgreader->find_by_path("CELLS/FIELD[@name='"+name+"']", true);
	return read_field<A>(curnode, Nc);
}
template<class A>
vector<A> GridReader::read_field(HMXML::Reader& rd, int num){
	TFieldInfo info(rd);
	auto content = preader->read_num_content(rd, num);
	if (info.type == "int") return content.convert_data<int, A>();
	if (info.type == "float") return content.convert_data<float, A>();
	if (info.type == "double") return content.convert_data<double, A>();
	if (info.type == "char") return content.convert_data<char, A>();
	return vector<A>();
}

template<class A>
vector<vector<A>> GridReader::read_vertices_vecfield(std::string name){
	auto curnode = pgreader->find_by_path("VERTICES/FIELD[@name='"+name+"']", true);
	return read_vecfield<A>(curnode, Nv);
}
template<class A>
vector<vector<A>> GridReader::read_edges_vecfield(std::string name){
	auto curnode = pgreader->find_by_path("EDGES/FIELD[@name='"+name+"']", true);
	return read_vecfield<A>(curnode, Ne);
}
template<class A>
vector<vector<A>> GridReader::read_faces_vecfield(std::string name){
	auto curnode = pgreader->find_by_path("FACES/FIELD[@name='"+name+"']", true);
	return read_vecfield<A>(curnode, Nc);
}
template<class A>
vector<vector<A>> GridReader::read_cells_vecfield(std::string name){
	auto curnode = pgreader->find_by_path("CELLS/FIELD[@name='"+name+"']", true);
	return read_vecfield<A>(curnode, Nc);
}
template<class A>
vector<vector<A>> GridReader::read_vecfield(HMXML::Reader& rd, int num){
	TFieldInfo info(rd);
	auto content = preader->read_num_content(rd, num);
	if (info.type == "int") return content.convert_vdata<int, A>();
	if (info.type == "float") return content.convert_vdata<float, A>();
	if (info.type == "double") return content.convert_vdata<double, A>();
	if (info.type == "char") return content.convert_vdata<char, A>();
	return vector<vector<A>>();
}




struct TReadHMG: public HMCallback::ExecutorBase{
	HMCB_SET_PROCNAME("Importing 3d grid from hmg");
	HMCB_SET_DEFAULT_DURATION(100);

	std::unique_ptr<GridReader> _run(HMXML::ReaderA* reader, HMXML::Reader* subnode);
};
extern HMCallback::FunctionWithCallback<TReadHMG> ReadHMG;

}}



#endif
