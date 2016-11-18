#include "hmg_import_grid3d.hpp"
using namespace HMGrid3D;
using namespace HMXML;

HMCallback::FunctionWithCallback<Import::TReadHMG> Import::ReadHMG;

std::unique_ptr<Import::GridReader> Import::TReadHMG::_run(HMXML::ReaderA* reader, HMXML::Reader* subnode){
	callback->silent_step_after(20, "Reading data", 4);
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
	size_t toteddim=0;
	for (auto& it: faceedge) toteddim+=it.size();

	//read faces->cells
	callback->subprocess_step_after(1);
	tmp = ret->pgreader->find_by_path("FACES/CELL_CONNECT", true);
	vector<int> facecell = ret->preader->read_num_content(tmp, ret->Nf*2).vec<int>();
	callback->subprocess_fin();

	callback->step_after(20, "Building serial grid");
	//constructing serialized grid
	ret->result.reset(new SGrid());
	std::swap(vert, ret->result->vert);
	std::swap(edgevert, ret->result->edges);
	ret->result->faces.resize(3*ret->Nf+toteddim);
	auto it=ret->result->faces.begin();
	for (size_t i=0; i<ret->Nf; ++i){
		size_t dim = faceedge[i].size();
		*it++ = dim;
		std::copy(faceedge[i].begin(), faceedge[i].end(), it);
		it += dim;
		*it++ = facecell[2*i];
		*it++ = facecell[2*i+1];
	}
	
	//boundary conditions
	try{
		vector<int> bt = ret->read_faces_field<int>("__boundary_types__");
		std::map<int, std::vector<int>> btmap;
		for (size_t i=0; i<bt.size(); ++i) if (bt[i]!=0){
			auto fnd = btmap.find(bt[i]);
			if (fnd == btmap.end()) fnd = btmap.emplace(bt[i], vector<int>()).first;
			fnd->second.push_back(i);
		}
		for (auto m: btmap) if (m.first != 0){
			ret->result->bnd.push_back(m.first);
			ret->result->bnd.push_back(m.second.size());
			std::copy(m.second.begin(), m.second.end(),
					std::back_inserter(ret->result->bnd));
		}
	} catch (const XmlElementNotFound&){}

	//additional serial information
	callback->step_after(20, "Additional data");
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
