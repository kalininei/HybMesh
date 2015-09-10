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
