#include "hmg_import_grid3d.hpp"
using namespace HMGrid3D;
using namespace HMXML;

HMCallback::FunctionWithCallback<Import::TReadHMG> Import::ReadHMG;

std::unique_ptr<Import::GridReader> Import::TReadHMG::_run(HMXML::ReaderA* reader, HMXML::Reader* subnode){
	callback->silent_step_after(30, "Reading data", 4);
	std::unique_ptr<Import::GridReader> ret(new GridReader(reader, subnode));
	//read dimensions
	ret->pgreader->value_int("N_VERTICES", ret->Nv, true);
	ret->pgreader->value_int("N_EDGES", ret->Ne, true);
	ret->pgreader->value_int("N_CELLS", ret->Nc, true);
	ret->pgreader->value_int("N_FACES", ret->Nf, true);

	//read vertices
	callback->subprocess_step_after(1);
	Reader tmp = ret->pgreader->find_by_path("VERTICES/COORDS", true);
	vector<double> vert = ret->preader->read_num_content(tmp, 3*ret->Nv).vec<double>();

	//read edges->vertices
	callback->subprocess_step_after(1);
	tmp = ret->pgreader->find_by_path("EDGES/VERT_CONNECT", true);
	vector<int> edgevert = ret->preader->read_num_content(tmp, 2*ret->Ne).vec<int>();

	//read faces->edges
	callback->subprocess_step_after(1);
	tmp = ret->pgreader->find_by_path("FACES/EDGE_CONNECT", true);
	vector<vector<int>> faceedge = ret->preader->read_num_content(tmp, ret->Nf).vecvec<int>();

	//read faces->cells
	callback->subprocess_step_after(1);
	tmp = ret->pgreader->find_by_path("FACES/CELL_CONNECT", true);
	vector<int> facecell = ret->preader->read_num_content(tmp, ret->Nf*2).vec<int>();
	callback->subprocess_fin();

	callback->step_after(30, "Building serial grid");
	//constructing serialized grid
	ret->result.reset(new SGrid());
	std::swap(vert, ret->result->vert);
	std::swap(edgevert, ret->result->edges);
	std::swap(faceedge, ret->result->face_edge);
	std::swap(facecell, ret->result->face_cell);
	//boundary conditions
	try{
		ret->result->btypes = ret->read_faces_field<int>("__boundary_types__");
	} catch (const XmlElementNotFound&){
		ret->result->btypes.resize(ret->Nf, 0);
	}
	//additional serial information
	ret->result->supplement();

	//unserial
	callback->step_after(40, "Assembling grid");
	ret->result->actualize_data();
	
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



// ======================= Surface reader
HMCallback::FunctionWithCallback<Import::TReadHMC> Import::ReadHMC;

std::unique_ptr<Import::SurfaceReader> Import::TReadHMC::_run(HMXML::ReaderA* reader, HMXML::Reader* subnode){
	callback->silent_step_after(30, "Reading data", 3);
	std::unique_ptr<Import::SurfaceReader> ret(new SurfaceReader(reader, subnode));
	//read dimensions
	ret->psreader->value_int("N_VERTICES", ret->Nv, true);
	ret->psreader->value_int("N_EDGES", ret->Ne, true);
	ret->psreader->value_int("N_FACES", ret->Nf, true);

	//read vertices
	callback->subprocess_step_after(1);
	Reader tmp = ret->psreader->find_by_path("VERTICES/COORDS", true);
	vector<double> vert = ret->preader->read_num_content(tmp, 3*ret->Nv).vec<double>();

	//read edges->vertices
	callback->subprocess_step_after(1);
	tmp = ret->psreader->find_by_path("EDGES/VERT_CONNECT", true);
	vector<int> edgevert = ret->preader->read_num_content(tmp, 2*ret->Ne).vec<int>();

	//read faces->edges
	callback->subprocess_step_after(1);
	tmp = ret->psreader->find_by_path("FACES/EDGE_CONNECT", true);
	vector<vector<int>> faceedge = ret->preader->read_num_content(tmp, ret->Nf).vecvec<int>();

	callback->step_after(30, "Building serial surface");
	//constructing serialized grid
	ret->result.reset(new SSurface());
	std::swap(vert, ret->result->vert);
	std::swap(edgevert, ret->result->edges);
	std::swap(faceedge, ret->result->face_edge);
	//boundary conditions
	try{
		ret->result->btypes = ret->read_faces_field<int>("__boundary_types__");
	} catch (const XmlElementNotFound&){
		ret->result->btypes.resize(ret->Nf, 0);
	}
	//additional serial information
	ret->result->supplement();

	//unserial
	callback->step_after(40, "Assembling surface");
	ret->result->actualize_data();
	
	return ret;
}

Import::SurfaceReader::TFieldInfo::TFieldInfo(HMXML::Reader& field){
	name = field.attribute(".", "name");
	type = field.attribute(".", "type");
	std::string sdim = field.attribute(".", "dim");

	if (name=="" || type=="") throw std::runtime_error(
			"not enough attributes in FIELD node");
	dim = 1;
	if (sdim=="variable") dim=-1;
	if (sdim!="") dim=atoi(sdim.c_str());
}
