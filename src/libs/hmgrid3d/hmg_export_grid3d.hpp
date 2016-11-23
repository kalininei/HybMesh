#ifndef HMG_EXPORT_GRID3D_HPP
#define HMG_EXPORT_GRID3D_HPP
#include "serialize_grid3d.hpp"
#include "hmxmlreader.hpp"

namespace HMGrid3D{ namespace Export{

struct GridWriter{
	//holds data

	GridWriter(const SGrid& g,
		HMXML::ReaderA* writer,
		HMXML::Reader* subnode,
		std::string gridname, std::string tp);

	template<class A>
	void AddVertexData(std::string fieldname, const A& data, bool binary);
	template<class A>
	void AddEdgeData(std::string fieldname, const A& data, bool binary);
	template<class A>
	void AddFaceData(std::string fieldname, const A& data, bool binary);
	template<class A>
	void AddCellData(std::string fieldname, const A& data, bool binary);

	void AddFaceVertexConnectivity();
	void AddCellFaceConnectivity();
	void AddCellVertexConnectivity();
	void AddLinFemConnectivity();

	template<class A> bool is_binary(){return false;}
private:
	virtual void data_changed(){}
	template<class A>
	static void add_field(HMXML::Reader& subnode, std::string fieldname, const vector<A>& data,
			bool binary, HMXML::ReaderA& writer);

	HMXML::ReaderA* pwriter;
	const SGrid* grid;
	HMXML::Reader gwriter, vwriter, ewriter, fwriter, cwriter;
	std::string __tp;
};

template<> inline bool GridWriter::is_binary<int>(){return __tp == "bin";}
template<> inline bool GridWriter::is_binary<char>(){return __tp == "bin";}
template<> inline bool GridWriter::is_binary<double>(){return __tp != "ascii";}
template<> inline bool GridWriter::is_binary<float>(){return __tp != "ascii";}

template<class A>
void GridWriter::add_field(HMXML::Reader& subnode, std::string fieldname, const vector<A>& data,
		bool binary, HMXML::ReaderA& writer){
	HMXML::Reader wr = subnode.find_by_path("FIELD[@name='"+fieldname+"']");
	if (wr) wr.unlink_node();
	wr = subnode.new_child("FIELD");
	wr.new_attribute("name", fieldname);
	writer.set_num_content(data, wr, binary);
}

template<class A>
void GridWriter::AddVertexData(std::string fieldname, const A& data, bool binary){
	data_changed();
	GridWriter::add_field<typename A::value_type>(vwriter, fieldname, data, binary, *pwriter);
}
template<class A>
void GridWriter::AddEdgeData(std::string fieldname, const A& data, bool binary){
	data_changed();
	GridWriter::add_field<typename A::value_type>(ewriter, fieldname, data, binary, *pwriter);
}
template<class A>
void GridWriter::AddCellData(std::string fieldname, const A& data, bool binary){
	data_changed();
	GridWriter::add_field<typename A::value_type>(cwriter, fieldname, data, binary, *pwriter);
}
template<class A>
void GridWriter::AddFaceData(std::string fieldname, const A& data, bool binary){
	data_changed();
	GridWriter::add_field<typename A::value_type>(fwriter, fieldname, data, binary, *pwriter);
}

}}



#endif
