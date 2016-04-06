#include "serialize_grid3d.hpp"
#include "hmcallback.hpp"
#include "addalgo.hpp"

using namespace HMGrid3D;

//static functors declaration
HMCallback::FunctionWithCallback<SimpleSerialize::ConvertExe>
SimpleSerialize::Convert;
HMCallback::FunctionWithCallback<ExtendedSimpleSerialize::ConvertExe>
ExtendedSimpleSerialize::Convert;

auto SimpleSerialize::ConvertExe::_run(const Grid& g)->TRet1{
	Grid::Talldata ad;
	return _run(g, ad);
}
auto SimpleSerialize::ConvertExe::_run(const Grid& g, Grid::Talldata& _alldata)->TRet1{
	SimpleSerialize ret;

	callback.step_after(40, "Gathering data");
	_alldata = g.alldata();
	auto& pvert = std::get<0>(_alldata);
	auto& pedges = std::get<1>(_alldata);
	auto& pfaces = std::get<2>(_alldata);
	auto& pcells = std::get<3>(_alldata);

	//find vertices of non-linear edges and remove'em from vert
	callback.step_after(5, "Simplify data");
	std::set<shared_ptr<Vertex>> not_used_vert;
	for (auto& e: pedges) if (e->vertices.size()>2){
		not_used_vert.insert(e->vertices.begin()+1, e->vertices.end()-1);
	}
	if (not_used_vert.size()>0){
		std::set<int> not_used_vert_ind;
		auto b = pvert.begin();
		for (auto v: not_used_vert){
			auto fnd = std::find(b, pvert.end(), v);
			assert(fnd != pvert.end());
			not_used_vert_ind.insert(fnd - pvert.begin());
			b = fnd;
		}
		aa::remove_entries(pvert, not_used_vert_ind);
	}

	//write primitives counts
	ret.n_vert = pvert.size();
	ret.n_edges = pedges.size();
	ret.n_faces = pfaces.size();
	ret.n_cells = pcells.size();

	callback.step_after(5, "Serialize points");
	ret.vert.reserve(3*ret.n_vert);
	for (auto& v: pvert){
		ret.vert.push_back(v->x);
		ret.vert.push_back(v->y);
		ret.vert.push_back(v->z);
	}

	callback.silent_step_after(10, "Serialize edges", 1);
	ret.edges.reserve(2*ret.n_edges);
	auto v_indexer = aa::shp_container_indexer(pvert);
	callback.subprocess_step_after(0.5);
	v_indexer.convert();
	callback.subprocess_step_after(0.1);
	for (auto& e: pedges){
		ret.edges.push_back(v_indexer.index(e->first()));
		ret.edges.push_back(v_indexer.index(e->last()));
	}
	callback.subprocess_step_after(0.4);
	v_indexer.restore();
	callback.subprocess_fin();

	callback.silent_step_after(15, "Serialize cells", 1);
	int szc = pcells.size();
	for (auto& c: pcells) szc += c->n_faces();
	ret.cells.reserve(szc);
	auto f_indexer = aa::shp_container_indexer(pfaces);
	callback.subprocess_step_after(0.5);
	f_indexer.convert();
	callback.subprocess_step_after(0.1);
	for (auto& c: pcells){
		ret.cells.push_back(c->n_faces());
		for (auto& f: c->faces){
			ret.cells.push_back(f_indexer.index(f));
		}
	}
	callback.subprocess_step_after(0.4);
	f_indexer.restore();
	callback.subprocess_fin();

	callback.silent_step_after(20, "Serialize faces", 1);
	int szf = pfaces.size()*3;
	for (auto& f: pfaces) szf+=f->n_edges();
	ret.faces.reserve(szf);
	auto e_indexer = aa::shp_container_indexer(pedges);
	auto c_indexer = aa::shp_container_indexer(pcells);
	callback.subprocess_step_after(0.25);
	e_indexer.convert();
	callback.subprocess_step_after(0.25);
	c_indexer.convert();
	callback.subprocess_step_after(0.1);
	for (auto& f: pfaces){
		ret.faces.push_back(f->n_edges());
		for (auto& e: f->edges) ret.faces.push_back(e_indexer.index(e));
		int cl = (f->left) ? c_indexer.index(f->left) : -1;
		int cr = (f->right) ? c_indexer.index(f->right) : -1;
		ret.faces.push_back(cl);
		ret.faces.push_back(cr);
	}
	callback.subprocess_step_after(0.2);
	e_indexer.restore();
	callback.subprocess_step_after(0.2);
	c_indexer.restore();
	callback.subprocess_fin();

	//write boundary conditions
	callback.step_after(5, "Serialize boundary");
	std::map<int, vector<int>> bc;
	int i = 0;
	for (auto& f: pfaces) {
		if (f->is_boundary() && f->boundary_type != 0){
			auto fnd = bc.find(f->boundary_type);
			if (fnd == bc.end()){
				auto emp = bc.emplace(f->boundary_type, std::vector<int>());
				fnd = emp.first;
			}
			fnd->second.push_back(i);
		}
		++i;
	}
	int szbnd = bc.size()*2;
	for (auto& b: bc) szbnd += b.second.size();
	ret.bnd.reserve(szbnd);
	for (auto& b: bc){
		ret.bnd.push_back(b.first);
		ret.bnd.push_back(b.second.size());
		ret.bnd.insert(ret.bnd.end(), b.second.begin(), b.second.end());
	}
	return ret;
}

//ExtendedSimpleSerialize::Convert
auto ExtendedSimpleSerialize::ConvertExe::_run(const Grid& g)->TRet1{
	ExtendedSimpleSerialize ret;

	auto r1 = SimpleSerialize::Convert.MoveCallback(callback,
			g, ret._alldata);
	ret.SimpleSerialize::operator=(std::move(r1));
	
	callback.step_after(10, "Extended data assembling");
	int k;
	ret.icell.reserve(ret.n_cells+1);
	k = 0;
	for (int i=0; i<ret.n_cells; ++i){
		ret.icell.push_back(k);
		k += ret.cells[k]+1;
	}
	ret.icell.push_back(k);
	ret.iface.reserve(ret.n_faces+1);
	k = 0;
	for (int i=0; i<ret.n_faces; ++i){
		ret.iface.push_back(k);
		k += ret.faces[k] + 3;
	}
	ret.iface.push_back(k);

	return ret;
}
