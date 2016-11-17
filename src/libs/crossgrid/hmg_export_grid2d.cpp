#include "hmg_export_grid2d.hpp"
#include <sstream>
#include <fstream>

using namespace GGeom;
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

Export::GridWriter::GridWriter(const GridGeom& g,
		HMXML::ReaderA* writer,
		HMXML::Reader* subnode,
		std::string gridname,
		std::string tp): pwriter(writer), grid(&g){
	__tp = tp;
	
	//supplementary data
	std::set<Edge> sedges = grid->get_edges();
	vector<Edge> edges(sedges.begin(), sedges.end());
	vector<double> pcoords(2*grid->n_points());
	for (int i=0; i<grid->n_points(); ++i){
		pcoords[2*i]   = grid->get_point(i)->x;
		pcoords[2*i+1] = grid->get_point(i)->y;
	}
	vector<int> edgeconnect(edges.size()*2), edgecellconnect(edges.size()*2);
	for (int i=0; i<edges.size(); ++i){
		edgeconnect[2*i] = edges[i].p1;
		edgeconnect[2*i+1] = edges[i].p2;
		edgecellconnect[2*i] = edges[i].cell_left;
		edgecellconnect[2*i+1] = edges[i].cell_right;
	}

	//create xml structure
	gwriter = subnode->new_child("GRID2D");
	gwriter.new_attribute("name", gridname);

	//inforamtion
	gwriter.new_child("N_VERTICES").set_content(std::to_string(grid->n_points()));
	gwriter.new_child("N_EDGES").set_content(std::to_string(edges.size()));
	gwriter.new_child("N_CELLS").set_content(std::to_string(grid->n_cells()));

	//vertices
	vwriter = gwriter.new_child("VERTICES");
	auto coordswriter = vwriter.new_child("COORDS");
	writer->set_num_content(pcoords, coordswriter, __tp!="ascii");

	//edges
	ewriter = gwriter.new_child("EDGES");
	auto vconnectwriter = ewriter.new_child("VERT_CONNECT");
	writer->set_num_content(edgeconnect, vconnectwriter, __tp=="bin");

	auto cconnectwriter = ewriter.new_child("CELL_CONNECT");
	writer->set_num_content(edgecellconnect, cconnectwriter, __tp=="bin");

	//cells
	cwriter = gwriter.new_child("CELLS");
}

void Export::GridWriter::AddCellVertexConnectivity(){
	data_changed();
	vector<vector<int>> data(grid->n_cells());
	for (int i=0; i<grid->n_cells(); ++i){
		auto cell = grid->get_cell(i);
		data[i].resize(cell->dim());
		for (int j=0; j<cell->dim(); ++j){
			data[i][j] = cell->get_point(j)->get_ind();
		}
	}
	AddCellData("__cell_vertices__", data, __tp=="bin");
}
void Export::GridWriter::AddCellEdgeConnectivity(){
	data_changed();
	vector<vector<int>> data(grid->n_cells());
	std::set<Edge> sedges = grid->get_edges();
	std::map<Edge, int> edges;
	auto it = sedges.begin();
	for (int i=0; i<sedges.size(); ++i) edges.emplace(*it++, i);
	for (int i=0; i<grid->n_cells(); ++i){
		const Cell* cell = grid->get_cell(i);
		data[i].resize(cell->dim());
		for (int j=0; j<cell->dim(); ++j){
			int p1 = cell->get_point(j)->get_ind();
			int p2 = cell->get_point(j+1)->get_ind();
			auto fnd = edges.find(Edge(p1, p2));
			assert(fnd != edges.end());
			data[i][j] = fnd->second;
		}
	}
	AddCellData("__cell_edges__", data, __tp=="bin");
}

void Export::GridWriter::AddVertexData(std::string fieldname, const vector<double>& data, bool binary){
	data_changed();
	add_field<double>(vwriter, fieldname, data, binary, *pwriter);
}
void Export::GridWriter::AddVertexData(std::string fieldname, const vector<float>& data, bool binary){
	data_changed();
	add_field<float>(vwriter, fieldname, data, binary, *pwriter);
}
void Export::GridWriter::AddVertexData(std::string fieldname, const vector<char>& data, bool binary){
	data_changed();
	add_field<char>(vwriter, fieldname, data, binary, *pwriter);
}
void Export::GridWriter::AddVertexData(std::string fieldname, const vector<int>& data, bool binary){
	data_changed();
	add_field<int>(vwriter, fieldname, data, binary, *pwriter);
}

void Export::GridWriter::AddEdgeData(std::string fieldname, const vector<double>& data, bool binary){
	data_changed();
	add_field<double>(ewriter, fieldname, data, binary, *pwriter);
}
void Export::GridWriter::AddEdgeData(std::string fieldname, const vector<float>& data, bool binary){
	data_changed();
	add_field<float>(ewriter, fieldname, data, binary, *pwriter);
}
void Export::GridWriter::AddEdgeData(std::string fieldname, const vector<char>& data, bool binary){
	data_changed();
	add_field<char>(ewriter, fieldname, data, binary, *pwriter);
}
void Export::GridWriter::AddEdgeData(std::string fieldname, const vector<int>& data, bool binary){
	data_changed();
	add_field<int>(ewriter, fieldname, data, binary, *pwriter);
}

void Export::GridWriter::AddCellData(std::string fieldname, const vector<double>& data, bool binary){
	data_changed();
	add_field<double>(cwriter, fieldname, data, binary, *pwriter);
}
void Export::GridWriter::AddCellData(std::string fieldname, const vector<float>& data, bool binary){
	data_changed();
	add_field<float>(cwriter, fieldname, data, binary, *pwriter);
}
void Export::GridWriter::AddCellData(std::string fieldname, const vector<char>& data, bool binary){
	data_changed();
	add_field<char>(cwriter, fieldname, data, binary, *pwriter);
}
void Export::GridWriter::AddCellData(std::string fieldname, const vector<int>& data, bool binary){
	data_changed();
	add_field<int>(cwriter, fieldname, data, binary, *pwriter);
}
void Export::GridWriter::AddCellData(std::string fieldname, const vector<vector<double>>& data, bool binary){
	data_changed();
	add_field<vector<double>>(cwriter, fieldname, data, binary, *pwriter);
}
void Export::GridWriter::AddCellData(std::string fieldname, const vector<vector<float>>& data, bool binary){
	data_changed();
	add_field<vector<float>>(cwriter, fieldname, data, binary, *pwriter);
}
void Export::GridWriter::AddCellData(std::string fieldname, const vector<vector<char>>& data, bool binary){
	data_changed();
	add_field<vector<char>>(cwriter, fieldname, data, binary, *pwriter);
}
void Export::GridWriter::AddCellData(std::string fieldname, const vector<vector<int>>& data, bool binary){
	data_changed();
	add_field<vector<int>>(cwriter, fieldname, data, binary, *pwriter);
}
void Export::GridWriter::AddVertexData(std::string fieldname, const vector<vector<double>>& data, bool binary){
	data_changed();
	add_field<vector<double>>(vwriter, fieldname, data, binary, *pwriter);
}
void Export::GridWriter::AddVertexData(std::string fieldname, const vector<vector<float>>& data, bool binary){
	data_changed();
	add_field<vector<float>>(vwriter, fieldname, data, binary, *pwriter);
}
void Export::GridWriter::AddVertexData(std::string fieldname, const vector<vector<char>>& data, bool binary){
	data_changed();
	add_field<vector<char>>(vwriter, fieldname, data, binary, *pwriter);
}
void Export::GridWriter::AddVertexData(std::string fieldname, const vector<vector<int>>& data, bool binary){
	data_changed();
	add_field<vector<int>>(vwriter, fieldname, data, binary, *pwriter);
}

void Export::GridWriter::AddEdgeData(std::string fieldname, const vector<vector<double>>& data, bool binary){
	data_changed();
	add_field<vector<double>>(ewriter, fieldname, data, binary, *pwriter);
}
void Export::GridWriter::AddEdgeData(std::string fieldname, const vector<vector<float>>& data, bool binary){
	data_changed();
	add_field<vector<float>>(ewriter, fieldname, data, binary, *pwriter);
}
void Export::GridWriter::AddEdgeData(std::string fieldname, const vector<vector<char>>& data, bool binary){
	data_changed();
	add_field<vector<char>>(ewriter, fieldname, data, binary, *pwriter);
}
void Export::GridWriter::AddEdgeData(std::string fieldname, const vector<vector<int>>& data, bool binary){
	data_changed();
	add_field<vector<int>>(ewriter, fieldname, data, binary, *pwriter);
}


// ==================================== GridWriter
Export::GridHMG::GridHMG(const GridGeom& g, std::string fn,
		std::string gridname, std::string tp):
		WOwner(),
		GridWriter(g, writer.get(), writer.get(), gridname, tp),
		__filename(fn), __flushed(false){}
		
Export::GridHMG::GridHMG(const GridGeom& g, std::string fn,
		shared_ptr<HMXML::ReaderA> writer,
		std::string gridname, std::string tp):
		WOwner(writer),
		GridWriter(g, writer.get(), writer.get(), gridname, tp),
		__filename(fn), __flushed(false){}

Export::GridHMG::~GridHMG(){
	Flush();
	if (writer.unique()) writer->Free();
}

void Export::GridHMG::Flush(){
	if (!__flushed && __filename != ""){
		__flushed = true;
		writer->write(__filename);
	}
}

// ==================================== MultipleGrids
Export::MultipleGridsHMG::MultipleGridsHMG(
		vector<const GridGeom*> grids,
		vector<std::string> names,
		std::string fn, std::string tp){
	if (grids.size() != names.size()){
		throw std::runtime_error("Number of grids doesn't equal number of names");
	}
	int N=grids.size();
	subs.reserve(N);
	writer.reset(new ReaderA(ReaderA::create("HybMeshData")));
	for (int i=0; i<grids.size(); ++i){
		subs.push_back(GridHMG(*grids[i], "", writer, names[i], tp));
	}
	subs[0].__filename=fn;
}
Export::MultipleGridsHMG::~MultipleGridsHMG(){
	Flush();
	if (writer.unique()) writer->Free();
}
void Export::MultipleGridsHMG::Flush(){
	if (!__flushed()){
		subs[0].Flush();
		for (int i=0; i<subs.size(); ++i) subs[i].__flushed=true;
	}
}
bool Export::MultipleGridsHMG::__flushed(){
	for (int i=0; i<subs.size(); ++i){
		if (subs[i].__flushed == false) return false;
	}
	return true;
}
