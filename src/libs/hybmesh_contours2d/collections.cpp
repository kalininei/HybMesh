#include "collections.hpp"
#include <fstream>
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
	std::set<Point*> tmp;
	std::for_each(data.begin(), data.end(), [&tmp](const Tentry& e){
			tmp.insert(e->pstart);
			tmp.insert(e->pend);
	});
	return vector<Point*>(tmp.begin(), tmp.end());
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
		if (dnew < dist){
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
	return BoundingBox::Build(ap.begin(), ap.end());
}



