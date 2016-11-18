#include "hmg_export_grid3d.hpp"

using namespace HMGrid3D;

Export::GridWriter::GridWriter(const SGrid& g,
		HMXML::ReaderA* writer,
		HMXML::Reader* subnode,
		std::string gridname, std::string tp){
	__tp = tp;
	grid=&g;
	pwriter=writer;

	//supplementary data
	vector<vector<int>> face_edge(g.n_faces);
	vector<int> face_cell(g.n_faces*2);
	auto it = g.faces.begin();
	auto it2 = face_cell.begin();
	for (int i=0; i<g.n_faces; ++i){
		int d = *it++;
		face_edge[i].resize(d);
		std::copy(it, it+d, face_edge[i].begin());
		it+=d;
		*it2++ = *it++;
		*it2++ = *it++;
	}
	
	//create xml structure
	gwriter = subnode->new_child("GRID3D");
	gwriter.new_attribute("name", gridname);

	//inforamtion
	gwriter.new_child("N_VERTICES").set_content(std::to_string(grid->n_vert));
	gwriter.new_child("N_FACES").set_content(std::to_string(grid->n_faces));
	gwriter.new_child("N_EDGES").set_content(std::to_string(grid->n_edges));
	gwriter.new_child("N_CELLS").set_content(std::to_string(grid->n_cells));

	//vertices
	vwriter = gwriter.new_child("VERTICES");
	auto coordswriter = vwriter.new_child("COORDS");
	writer->set_num_content(grid->vert, coordswriter, is_binary<double>());

	//edges
	ewriter = gwriter.new_child("EDGES");
	auto vconnectwriter = ewriter.new_child("VERT_CONNECT");
	writer->set_num_content(grid->edges, vconnectwriter, is_binary<int>());

	//faces
	fwriter = gwriter.new_child("FACES");
	auto econnectwriter = fwriter.new_child("EDGE_CONNECT");
	writer->set_num_content(face_edge, econnectwriter, is_binary<int>());

	auto cconnectwriter = fwriter.new_child("CELL_CONNECT");
	writer->set_num_content(face_cell, cconnectwriter, is_binary<int>());

	//cells
	cwriter = gwriter.new_child("CELLS");
	
	//boundary types
	std::vector<int> bt(g.n_faces); 
	for (size_t i=0; i<bt.size(); ++i) bt[i] = g.vfaces[i]->boundary_type;
	int minv = *std::min_element(bt.begin(), bt.end());
	int maxv = *std::max_element(bt.begin(), bt.end());
	if (minv == maxv && minv == 0) return;
	if (minv >-128 && maxv < 128){
		std::vector<char> btchar(bt.begin(), bt.end());
		AddFaceData("__boundary_types__", btchar, is_binary<char>());
	} else {
		AddFaceData("__boundary_types__", bt, is_binary<int>());
	}
}

void Export::GridWriter::AddFaceVertexConnectivity(){
	enumerate_ids_pvec(grid->vvert);
	std::vector<std::vector<int>> face_vertex(grid->n_faces);
	for (size_t i=0; i<face_vertex.size(); ++i){
		auto fv = grid->vfaces[i]->sorted_vertices();
		face_vertex[i].resize(fv.size());
		for (size_t j=0; j<fv.size(); ++j){
			face_vertex[i][j] = fv[j]->id;
		}
	}
	AddFaceData("__face_vertices__", face_vertex, is_binary<int>());
}
void Export::GridWriter::AddCellFaceConnectivity(){
	enumerate_ids_pvec(grid->vfaces);
	std::vector<std::vector<int>> cell_face(grid->n_cells);
	for (size_t i=0; i<cell_face.size(); ++i){
		auto& cf = grid->vcells[i]->faces;
		cell_face[i].resize(cf.size());
		for (size_t j=0; j<cf.size(); ++j){
			cell_face[i][j] = cf[j]->id;
		}
	}
	AddCellData("__cell_faces__", cell_face, is_binary<int>());
}
void Export::GridWriter::AddCellVertexConnectivity(){
	enumerate_ids_pvec(grid->vvert);
	std::vector<std::vector<int>> cell_vertex(grid->n_cells);
	std::vector<int> tmp;
	size_t i=0;
	for (auto& c: grid->vcells){
		tmp.resize(0);
		for (auto& f: c->faces)
		for (auto& e: f->edges)
		for (auto& v: e->vertices){
			tmp.push_back(v->id);
		}
		auto itback = std::unique(tmp.begin(), tmp.end());
		std::copy(tmp.begin(), itback, std::back_inserter(cell_vertex[i++]));
	}
	AddCellData("__cell_vertices__", cell_vertex, is_binary<int>());
}
