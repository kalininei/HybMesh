#include "collections.hpp"
#include <fstream>

using namespace HMCont2D;

void ECollection::SaveVtk(const ECollection& ec, const char* fn){
	std::ofstream fs(fn);
	fs<<"# vtk DataFile Version 3.0"<<std::endl;
	fs<<"EdgesData"<<std::endl;
	fs<<"ASCII"<<std::endl;
	//Points
	vector<Point*> pts = ec.all_points();
	fs<<"DATASET UNSTRUCTURED_GRID"<<std::endl;
	fs<<"POINTS "<<pts.size()<< " float"<<std::endl;
	for (int i=0; i<pts.size(); ++i){
		fs<<pts[i]->x<<" "<<pts[i]->y<<" 0"<<std::endl;
	}
	//Cells
	fs<<"CELLS  "<<ec.size()<<"   "<<3*ec.size()<<std::endl;
	for (int i=0; i<ec.size(); ++i){
		int i1 = std::find(pts.begin(), pts.end(), ec.pvalue(i)->pstart) - pts.begin();
		int i2 = std::find(pts.begin(), pts.end(), ec.pvalue(i)->pend) - pts.begin();
		fs<<2<<" "<<i1<<" "<<i2<<std::endl;
	}
	fs<<"CELL_TYPES  "<<ec.size()<<std::endl;
	for (int i=0;i<ec.size();++i) fs<<3<<std::endl;
	fs.close();
}
