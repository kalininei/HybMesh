#include "export2d_hm.hpp"
#include "contour.hpp"

using namespace HM2D;
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

Export::EColWriter::EColWriter(const EdgeData& c,
		HMXML::ReaderA* writer,
		HMXML::Reader* subnode,
		std::string contname,
		std::string tp): pwriter(writer), cont(&c){
	__tp = tp;
	
	//supplementary data
	auto ap = AllVertices(c);
	vector<double> pcoords(ap.size()*2);
	for (int i=0; i<ap.size(); ++i){
		pcoords[2*i] = ap[i]->x;
		pcoords[2*i+1] = ap[i]->y;
	}
	vector<int> edgeconnect; edgeconnect.reserve(c.size()*2);
	auto _indexer = aa::ptr_container_indexer(ap);
	_indexer.convert();
	for (auto e: *cont){
		edgeconnect.push_back(_indexer.index(e->first().get()));
		edgeconnect.push_back(_indexer.index(e->last().get()));
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

	//boundary types
	vector<int> bt(c.size(), 0);
	for (int i=0; i<c.size(); ++i) bt[i] = c[i]->boundary_type;
	//calculate min and max to select best storage type
	int minv = *std::min_element(bt.begin(), bt.end());
	int maxv = *std::max_element(bt.begin(), bt.end());
	if (minv == maxv && minv == 0){
		//pass
	} else if (minv >-128 && maxv < 128){
		std::vector<char> btchar(bt.begin(), bt.end());
		AddEdgeData("__boundary_types__", btchar, is_binary<char>());
	} else {
		AddEdgeData("__boundary_types__", bt, is_binary<int>());
	}
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

Export::GridWriter::GridWriter(const GridData& g,
		HMXML::ReaderA* writer,
		HMXML::Reader* subnode,
		std::string gridname,
		std::string tp): pwriter(writer), grid(&g){
	g.enumerate_all();
	__tp = tp;
	
	//supplementary data
	vector<double> pcoords(2*grid->vvert.size());
	for (int i=0; i<grid->vvert.size(); ++i){
		pcoords[2*i]   = grid->vvert[i]->x;
		pcoords[2*i+1] = grid->vvert[i]->y;
	}
	vector<int> edgeconnect(grid->vedges.size()*2), edgecellconnect(grid->vedges.size()*2);
	for (int i=0; i<grid->vedges.size(); ++i){
		auto e = grid->vedges[i];
		edgeconnect[2*i] = e->first()->id;
		edgeconnect[2*i+1] = e->last()->id;
		int cl = -1, cr = -1;
		if (e->has_left_cell()) cl = e->left.lock()->id;
		if (e->has_right_cell()) cr = e->right.lock()->id;
		edgecellconnect[2*i] = cl;
		edgecellconnect[2*i+1] = cr;
	}

	//create xml structure
	gwriter = subnode->new_child("GRID2D");
	gwriter.new_attribute("name", gridname);

	//inforamtion
	gwriter.new_child("N_VERTICES").set_content(std::to_string(grid->vvert.size()));
	gwriter.new_child("N_EDGES").set_content(std::to_string(grid->vedges.size()));
	gwriter.new_child("N_CELLS").set_content(std::to_string(grid->vcells.size()));

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

	//boundary types
	vector<int> bt(g.vedges.size(), 0);
	for (int i=0; i<g.vedges.size(); ++i) bt[i] = g.vedges[i]->boundary_type;
	//calculate min and max to select best storage type
	int minv = *std::min_element(bt.begin(), bt.end());
	int maxv = *std::max_element(bt.begin(), bt.end());
	if (minv == maxv && minv == 0){
		//pass
	} else if (minv >-128 && maxv < 128){
		std::vector<char> btchar(bt.begin(), bt.end());
		AddEdgeData("__boundary_types__", btchar, is_binary<char>());
	} else {
		AddEdgeData("__boundary_types__", bt, is_binary<int>());
	}
}

void Export::GridWriter::AddCellVertexConnectivity(){
	data_changed();
	vector<vector<int>> data(grid->vcells.size());
	aa::enumerate_ids_pvec(grid->vvert);
	for (int i=0; i<grid->vcells.size(); ++i){
		auto cell = grid->vcells[i];
		data[i].resize(cell->edges.size());
		auto op = Contour::OrderedPoints(cell->edges);
		for (int j=0; j<op.size()-1; ++j){
			data[i][j] = op[j]->id;
		}
	}
	AddCellData("__cell_vertices__", data, __tp=="bin");
}
void Export::GridWriter::AddCellEdgeConnectivity(){
	data_changed();
	vector<vector<int>> data(grid->vcells.size());
	aa::enumerate_ids_pvec(grid->vedges);
	for (int i=0; i<grid->vcells.size(); ++i){
		const Cell* cell = grid->vcells[i].get();
		data[i].resize(cell->edges.size());
		for (int j=0; j<cell->edges.size(); ++j){
			data[i][j] = cell->edges[j]->id;
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
Export::GridHMG::GridHMG(const GridData& g, std::string fn,
		std::string gridname, std::string tp):
		WOwner(),
		GridWriter(g, writer.get(), writer.get(), gridname, tp),
		__filename(fn), __flushed(false){}
		
Export::GridHMG::GridHMG(const GridData& g, std::string fn,
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
		vector<const GridData*> grids,
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
