#include <fstream>
#include "tecplot_export_grid3d.hpp"
#include "surface_grid3d.hpp"
#include "debug_grid3d.hpp"

namespace hme = HMGrid3D::Export;
HMCallback::FunctionWithCallback<hme::TGridTecplot> hme::GridTecplot;
HMCallback::FunctionWithCallback<hme::TBoundaryTecplot> hme::BoundaryTecplot;

namespace {
typedef aa::PtrContainerIndexer<const ShpVector<HMGrid3D::Vertex>> TVertIndexer; 

struct SurfSerial{
	SurfSerial(HMGrid3D::Surface& srf, TVertIndexer& vrt){
		n_faces = srf.faces.size();
		//edge->nodes
		auto ae = srf.alledges();
		n_edges = ae.size();
		edges.reserve(n_edges*2);
		for (auto e: ae){
			edges.push_back(vrt.index(e->first()));
			edges.push_back(vrt.index(e->last()));
		}
		//edge->faces
		vector<bool> isleft_face_edge;
		for (auto f: srf.faces){
			auto eprev = f->edges.back();
			for (int i=0; i<f->n_edges(); ++i){
				auto e = f->edges[i];
				bool isleft = (e->first() == eprev->first() || e->first() == eprev->last());
				if (!f->has_right_cell()) isleft = !isleft;
				isleft_face_edge.push_back(isleft);
				eprev = e;
			}
		}
		auto _eindexer = aa::ptr_container_indexer(ae);
		vector<HMGrid3D::Face*> fleft(n_edges, 0), fright(n_edges, 0);
		_eindexer.convert();
		auto it = isleft_face_edge.begin();
		for (auto f: srf.faces){
			for (auto e: f->edges){
				if (*it++) fleft[_eindexer.index(e)] = f.get();
				else fright[_eindexer.index(e)] = f.get();
			}
		}
		_eindexer.restore();
		//to integer indicies
		edge_adj.resize(n_edges*2, -1);
		auto _findexer = aa::ptr_container_indexer(srf.faces);
		_findexer.convert();
		for (int i=0; i<n_edges; ++i){
			if (fleft[i] != 0) edge_adj[2*i] = _findexer.index(fleft[i]);
			if (fright[i] != 0) edge_adj[2*i+1] = _findexer.index(fright[i]);
		}
		_findexer.restore();
	}
	int n_edges, n_faces;
	vector<int> edges; //start_node_index, end_node_index for each face
	vector<int> edge_adj;//left face, right face for each edge

	void serialize_vertices(const ShpVector<HMGrid3D::Vertex>& vert){
		n_vert = vert.size();
		vertices.resize(3*n_vert);
		auto it = vertices.begin();
		for (int i=0; i<n_vert; ++i){
			auto p = vert[i];
			*it++ = p->x;
			*it++ = p->y;
			*it++ = p->z;
		}
	}
	int n_vert;
	vector<double> vertices;
};

template<int Step, int From, int RowN, class Func, class V>
void write_row_n(std::ostream& str, Func&& fun, const vector<V>& vals){
	static_assert(RowN > 1, "RowN>1");
	auto it = vals.begin() + From;
	int len = vals.size()/Step;
	int lasti = len/RowN;
	int resj = len - lasti*RowN;

	for (int i=0; i<lasti; ++i){
		for (int j=0; j<RowN-1; ++j){
			str<<fun(*it)<<" "; it += Step;
		}
		str<<fun(*it)<<std::endl; it+=Step;
	}
	for (int j=0; j<resj-1; ++j){
		str<<fun(*it)<<" "; it += Step;
	}
	if (resj>0) str<<fun(*it)<<std::endl; it+=Step;
}

};

void hme::TGridTecplot::_run(const GridData& g, std::string fn, BFun bnd_names){
	SGrid sg(g);
	return _run(sg, fn, bnd_names);
}
void hme::TGridTecplot::_run(const SGrid& ser, std::string fn, BFun bnames){
	callback->step_after(30, "Assembling connectivity");
	//face->nodes connectivity
	vector<vector<int>> face_nodes = ser.face_vertex();
	//total face connectivity
	int totalfn = 0;
	for (int i=0; i<face_nodes.size(); ++i) totalfn+=face_nodes[i].size();
	//face adjacents
	vector<int> left_cells, right_cells;
	{
		left_cells.reserve(ser.n_faces); right_cells.reserve(ser.n_faces);
		for (int i=0; i<ser.n_faces; ++i){
			left_cells.push_back(ser.face_cell[2*i]);
			right_cells.push_back(ser.face_cell[2*i+1]);
		}
	}
	//assembling surfaces
	std::map<int, HMGrid3D::Surface> surfaces_geom; 
	{
		for (auto& f: ser.vfaces) if (f->is_boundary()){
			int bt = f->boundary_type;
			auto fnd = surfaces_geom.find(bt);
			if (fnd == surfaces_geom.end()){
				fnd = surfaces_geom.emplace(bt, HMGrid3D::Surface()).first;
			}
			fnd->second.faces.push_back(f);
		}
	}
	//serializing surfaces
	callback->silent_step_after(20, "Serialize surfaces", surfaces_geom.size());
	std::map<int, SurfSerial> surfaces; 
	{
		auto _indexer = aa::ptr_container_indexer(ser.vvert);
		_indexer.convert();
		for (auto& m: surfaces_geom){
			callback->subprocess_step_after(1);
			surfaces.emplace(m.first, SurfSerial(m.second, _indexer));
		}
		_indexer.restore();
	}

	// ====== write to file:
	callback->silent_step_after(30, "Write to file", 10);
	std::ofstream of(fn);
	of.precision(10);
	//====== main header
	of<<"TITLE=\"Tecplot Export\""<<std::endl;
	of<<"VARIABLES=\"X\" \"Y\" \"Z\""<<std::endl;
	of<<"ZONE T=\"Grid\""<<std::endl;
	of<<"Nodes="<<ser.n_vert<<std::endl;
	of<<"Faces="<<ser.n_faces<<std::endl;
	of<<"Elements="<<ser.n_cells<<std::endl;
	of<<"ZONETYPE=FEPOLYHEDRON"<<std::endl;
	of<<"DATAPACKING=BLOCK"<<std::endl;
	of<<"TotalNumFaceNodes="<<totalfn<<std::endl;
	of<<"NumConnectedBoundaryFaces=0, TotalNumBoundaryConnections=0"<<std::endl;
	//points
	callback->subprocess_step_after(3);
	write_row_n<3, 0, 20>(of, [](double v){ return v; }, ser.vert);
	write_row_n<3, 1, 20>(of, [](double v){ return v; }, ser.vert);
	write_row_n<3, 2, 20>(of, [](double v){ return v; }, ser.vert);
	//face dims
	callback->subprocess_step_after(1);
	//for (int i=0; i<face_nodes.size(); ++i) of<<face_nodes[i].size()<<std::endl;
	write_row_n<1, 0, 20>(of, [](const vector<int>& v){ return v.size(); }, face_nodes);
	//face->nodes
	callback->subprocess_step_after(3);
	for (int i=0; i<face_nodes.size(); ++i){
		for (int j=0; j<face_nodes[i].size(); ++j) of<<face_nodes[i][j]+1<<" ";
		of<<std::endl;
	}
	//face left/right cells
	callback->subprocess_step_after(1);
	write_row_n<1, 0, 20>(of, [](const int& v){ return v + 1; }, left_cells);
	write_row_n<1, 0, 20>(of, [](const int& v){ return v + 1; }, right_cells);

	//====== boundary surfaces
	callback->subprocess_step_after(3);
	for (auto& s: surfaces){
		std::string name = bnames(s.first); 
		of<<"ZONE T=\""<<name<<"\""<<std::endl;
		of<<"D=(1 2 3)"<<std::endl;
		of<<"Faces="<<s.second.n_edges<<std::endl;
		of<<"Elements="<<s.second.n_faces<<std::endl;
		of<<"ZONETYPE=FEPOLYGON"<<std::endl;
		of<<"DATAPACKING=BLOCK"<<std::endl;
		of<<"NumConnectedBoundaryFaces=0, TotalNumBoundaryConnections=0"<<std::endl;
		//edges
		write_row_n<1, 0, 20>(of, [](const int& v){ return v+1; }, s.second.edges);
		//adjacents left
		write_row_n<2, 0, 20>(of, [](const int& v){ return v+1; }, s.second.edge_adj);
		//adjacents right
		write_row_n<2, 1, 20>(of, [](const int& v){ return v+1; }, s.second.edge_adj);
	}
	
	of.close();
}


void hme::TBoundaryTecplot::_run(const SGrid& g, std::string fn, BFun bnames){
	callback->step_after(30, "Assembling Surfaces");
	ShpVector<Face> af = g.vfaces;
	//assembling surfaces
	std::map<int, HMGrid3D::Surface> surfaces_geom; 
	{
		for (auto& f: af) if (f->is_boundary()){
			int bt = f->boundary_type;
			auto fnd = surfaces_geom.find(bt);
			if (fnd == surfaces_geom.end()){
				fnd = surfaces_geom.emplace(bt, HMGrid3D::Surface()).first;
			}
			fnd->second.faces.push_back(f);
		}
	}
	//serializing surfaces
	callback->silent_step_after(20, "Serialize surfaces", surfaces_geom.size());
	std::map<int, SurfSerial> surfaces; 
	{
		for (auto& m: surfaces_geom){
			callback->subprocess_step_after(1);
			const ShpVector<Vertex> allvert = m.second.allvertices();
			auto _indexer = aa::ptr_container_indexer(allvert);
			_indexer.convert();
			auto emp = surfaces.emplace(m.first, SurfSerial(m.second, _indexer));
			_indexer.restore();
			emp.first->second.serialize_vertices(allvert);
		}
	}

	// ====== write to file:
	callback->silent_step_after(30, "Write to file", surfaces.size());
	std::ofstream of(fn);
	of.precision(10);
	//====== main header
	of<<"TITLE=\"Tecplot Export\""<<std::endl;
	of<<"VARIABLES=\"X\" \"Y\" \"Z\""<<std::endl;
	//====== boundary surfaces
	for (auto& s: surfaces){
		callback->subprocess_step_after(1);
		std::string name = bnames(s.first); 
		of<<"ZONE T=\""<<name<<"\""<<std::endl;
		of<<"Nodes="<<s.second.n_vert<<std::endl;
		of<<"Faces="<<s.second.n_edges<<std::endl;
		of<<"Elements="<<s.second.n_faces<<std::endl;
		of<<"ZONETYPE=FEPOLYGON"<<std::endl;
		of<<"DATAPACKING=BLOCK"<<std::endl;
		of<<"NumConnectedBoundaryFaces=0, TotalNumBoundaryConnections=0"<<std::endl;
		//points
		write_row_n<3, 0, 20>(of, [](double v){ return v; }, s.second.vertices);
		write_row_n<3, 1, 20>(of, [](double v){ return v; }, s.second.vertices);
		write_row_n<3, 2, 20>(of, [](double v){ return v; }, s.second.vertices);
		//edges
		write_row_n<1, 0, 20>(of, [](const int& v){ return v+1; }, s.second.edges);
		//adjacents left
		write_row_n<2, 0, 20>(of, [](const int& v){ return v+1; }, s.second.edge_adj);
		//adjacents right
		write_row_n<2, 1, 20>(of, [](const int& v){ return v+1; }, s.second.edge_adj);
	}
	
	of.close();
}
