#ifndef HYBMESH_HM2D_IMPORT_HM_HPP
#define HYBMESH_HM2D_IMPORT_HM_HPP

#include "hmxmlreader.hpp"
#include "primitives2d.hpp"

namespace HM2D{namespace Import{

// ============================== Contour
struct EColReader{
	EColReader(HMXML::ReaderA* preader, HMXML::Reader* subnode);
	std::unique_ptr<HM2D::EdgeData> result;

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

// ===================================== Grid
struct GridReader{
	GridReader(HMXML::ReaderA* preader, HMXML::Reader* subnode);
	std::unique_ptr<GridData> result;

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
	template<class A> vector<A> read_cells_field(std::string name);

	template<class A> vector<vector<A>> read_vertices_vecfield(std::string name);
	template<class A> vector<vector<A>> read_edges_vecfield(std::string name);
	template<class A> vector<vector<A>> read_cells_vecfield(std::string name);

	int Nc, Ne, Nv;
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

struct WOwner{
	HMXML::Reader greader;
	shared_ptr<HMXML::ReaderA> reader;
	WOwner(std::string fn, std::string gname);
	~WOwner(){ if (reader.unique()) reader->Free(); }
};

struct GridHMG: private WOwner, public GridReader{
	GridHMG(std::string fn, std::string gname=""): WOwner(fn, gname), GridReader(reader.get(), &greader){}
};

}}



#endif
