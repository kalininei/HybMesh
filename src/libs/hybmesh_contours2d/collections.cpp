#include <fstream>
#include <unordered_map>
#include "collections.hpp"
using namespace HMCont2D;

ScaleBase PCollection::Scale01(PCollection& pc){
	return ScaleBase::p_doscale(pc.data.begin(), pc.data.end());

}
void PCollection::Scale(PCollection& pc, const ScaleBase& sc){
	sc.p_scale(pc.data.begin(), pc.data.end());
}
void PCollection::Unscale(PCollection& pc, const ScaleBase& sc){
	sc.p_unscale(pc.data.begin(), pc.data.end());
}

ScaleBase ECollection::Scale01(ECollection& pc){
	auto pa = pc.all_points();
	return ScaleBase::p_doscale(pa.begin(), pa.end());

}
void ECollection::Scale(ECollection& pc, const ScaleBase& sc){
	auto pa = pc.all_points();
	sc.p_scale(pa.begin(), pa.end());
}
void ECollection::Unscale(ECollection& pc, const ScaleBase& sc){
	auto pa = pc.all_points();
	sc.p_unscale(pa.begin(), pa.end());
}


std::tuple<Point*, int, double>
PCollection::FindClosestNode(const PCollection& dt, const Point& p){
	std::tuple<Point*, int, double> ret;
	Point*& p2 = std::get<0>(ret);
	int& ind = std::get<1>(ret);
	double& m = std::get<2>(ret);
	p2=0; ind=-1; m=1e200;
	for (int i=0; i<dt.size(); ++i){
		double m2 = Point::meas(p, *dt.data[i]);
		if (m2<m){
			m = m2;
			p2 = dt.data[i].get();
			ind = i;
		}
	}
	return ret;
}

vector<Point*> ECollection::all_points() const{
	vector<Point*> ret;
	auto ins = aa::unique_container_inserter(ret);
	std::for_each(data.begin(), data.end(), [&ins](const Tentry& e){
			ins.insert(e->pstart);
			ins.insert(e->pend);
	});
	return ret;
}
vector<std::array<int, 2>> ECollection::tab_edges_points() const{
	vector<Point*> ap = all_points();
	auto& ae = data;
	vector<std::array<int, 2>> edges_points(ae.size());
	auto _indexer = aa::ptr_container_indexer(ap);
	_indexer.convert();
	for (int i=0; i<ae.size(); ++i){
		edges_points[i][0] = _indexer.index(ae[i]->pstart);
		edges_points[i][1] = _indexer.index(ae[i]->pend);
	}
	_indexer.restore();
	return edges_points;
}
vector<vector<int>> ECollection::tab_points_edges() const{
	auto ep = tab_edges_points();
	int maxp = -1;
	for (auto& it: ep) maxp = std::max(std::max(maxp, it[0]), it[1]);
	vector<vector<int>> ret(maxp+1);
	for (int i=0; i<ep.size(); ++i){
		ret[ep[i][0]].push_back(i);
		ret[ep[i][1]].push_back(i);
	}
	return ret;
}
vector<vector<int>> ECollection::tab_edges_edges() const{
	auto pe = tab_points_edges();
	vector<vector<int>> ret(size());
	for (auto& it: pe)
	for (int i=1; i<it.size(); ++i)
	for (int j=0; j<i; ++j){
		ret[it[i]].push_back(it[j]);
		ret[it[j]].push_back(it[i]);
	}

	return ret;
}

Point* ECollection::FindClosestNode(const ECollection& dt, const Point& p){
	assert(dt.size() > 0);
	Point* res = dt.edge(0)->pstart;
	double d = Point::meas(p, *res);
	for (auto& e: dt.data){
		double d1 = Point::meas(p, *e->pstart);
		if (d1<d){ d = d1; res = e->pstart; }
		double d2 = Point::meas(p, *e->pend);
		if (d2<d){ d = d2; res = e->pend; }
	}
	return res;
}

std::tuple<Edge*, double, double, int>
ECollection::FindClosestEdge(const ECollection& dt, const Point& p){
	std::tuple<Edge*, double, double, int> ret;
	Edge*& edge = std::get<0>(ret);
	double& dist = std::get<1>(ret);
	double& ksi = std::get<2>(ret);
	int& ind = std::get<3>(ret);
	edge = 0; dist = 1e99; ksi = 1e99; ind = -1;

	double k;
	for (int i=0; i<dt.size(); ++i){
		auto& e = dt.data[i];
		double dnew = Point::meas_section(p, *e->pstart, *e->pend, k);
		//if (dnew < dist){
		if (dnew - dist < geps*geps){
			edge = e.get();
			dist = dnew;
			ksi = k;
			ind = i;
			if (dist<geps*geps) break;
		}
	}

	dist = sqrt(dist);
	return ret;
}

Point ECollection::ClosestPoint(const ECollection& dt, const Point& p){
	auto ce = FindClosestEdge(dt, p);
	if (ISEQ(std::get<2>(ce), 0)){
		return *std::get<0>(ce)->pstart;
	} else if (ISEQ(std::get<2>(ce), 1)){
		return *std::get<0>(ce)->pend;
	} else {
		Point* p1 = std::get<0>(ce)->pstart;
		Point* p2 = std::get<0>(ce)->pend;
		return Point::Weigh(*p1, *p2, std::get<2>(ce));
	}
}

void ECollection::SaveVtk(const ECollection& dt, const char* fn){
	std::ofstream fs(fn);
	fs<<"# vtk DataFile Version 3.0"<<std::endl;
	fs<<"Contour 2D"<<std::endl;
	fs<<"ASCII"<<std::endl;
	//Points
	vector<Point*> pts = dt.all_points();
	fs<<"DATASET UNSTRUCTURED_GRID"<<std::endl;
	fs<<"POINTS "<<pts.size()<< " float"<<std::endl;
	for (int i=0; i<pts.size(); ++i){
		fs<<pts[i]->x<<" "<<pts[i]->y<<" 0"<<std::endl;
	}
	//Cells
	fs<<"CELLS  "<<dt.size()<<"   "<<3*dt.size()<<std::endl;
	for (int i=0; i<dt.size(); ++i){
		int i1 = std::find(pts.begin(), pts.end(), dt.edge(i)->pstart) - pts.begin();
		int i2 = std::find(pts.begin(), pts.end(), dt.edge(i)->pend) - pts.begin();
		fs<<2<<" "<<i1<<" "<<i2<<std::endl;
	}
	fs<<"CELL_TYPES  "<<dt.size()<<std::endl;
	for (int i=0;i<dt.size();++i) fs<<3<<std::endl;
	fs.close();
}

void PCollection::SaveVtk(const PCollection& dt, const char* fn){
	std::ofstream fs(fn);
	fs<<"# vtk DataFile Version 3.0"<<std::endl;
	fs<<"Contour 2D"<<std::endl;
	fs<<"ASCII"<<std::endl;
	//Points
	auto& d = dt.data;
	fs<<"DATASET UNSTRUCTURED_GRID"<<std::endl;
	fs<<"POINTS "<<dt.size()<< " float"<<std::endl;
	for (int i=0; i<dt.size(); ++i){
		fs<<d[i]->x<<" "<<d[i]->y<<" 0"<<std::endl;
	}
	//Cells
	fs<<"CELLS  "<<dt.size()<<"   "<<2*dt.size()<<std::endl;
	for (int i=0; i<dt.size(); ++i){
		fs<<1<<" "<<i<<std::endl;
	}
	fs<<"CELL_TYPES  "<<dt.size()<<std::endl;
	for (int i=0;i<dt.size();++i) fs<<1<<std::endl;
	fs.close();
}


vector<double> ECollection::ELengths(const ECollection& dt){
	vector<double> ret; ret.reserve(dt.size());
	for (auto it = dt.begin(); it != dt.end(); ++it) ret.push_back((*it)->length());
	return ret;
}


BoundingBox ECollection::BBox(const ECollection& p, double eps){
	auto ap = p.all_points();
	assert(ap.size() > 0);
	auto ret = BoundingBox::Build(ap.begin(), ap.end());
	ret.widen(eps);
	return ret;
}

void ECollection::ReallocatePoints(PCollection& pcol){
	std::unordered_map<Point*, shared_ptr<Point>> inserter;
	for (auto e: data){
		assert(e->pstart != NULL && e->pend != NULL);
		auto fnd1 = inserter.find(e->pstart);
		if (fnd1 == inserter.end()){
			auto emp = inserter.emplace(e->pstart,
				shared_ptr<Point>(new Point(*e->pstart)));
			fnd1 = emp.first;
			pcol.add_value(fnd1->second);
		}
		auto fnd2 = inserter.find(e->pend);
		if (fnd2 == inserter.end()){
			auto emp = inserter.emplace(e->pend,
				shared_ptr<Point>(new Point(*e->pend)));
			fnd2 = emp.first;
			pcol.add_value(fnd2->second);
		}

		e->pstart = fnd1->second.get();
		e->pend = fnd2->second.get();
	}
}
