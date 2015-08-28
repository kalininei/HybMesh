#ifndef CROSSGRID_FILEPROC_H
#define CROSSGRID_FILEPROC_H
#include <chrono>
#include "bgeom2d.h"
#include "grid.h"
#include "wireframegrid.h"

// ---------------- debugging procedures
void save_vtk(const GridGeom* g, const char* fn);
void save_vtk(const PContour* c, const char* fn);
void save_vtk(const vector<PContour>& c, const char* fn);
void save_vtk(const vector<PContour>& c, const vector<double>& data, const char* fn);
void save_vtk(const PtsGraph* g, const char* fn);

void save_vtk(const ContoursCollection& c, const char* fn);
void save_vtk(const ContoursCollection* c, const char* fn);
void save_vtk(const ContoursCollection& c, const vector<double>& pdata, const char* fn);
void save_vtk(const GridGeom& g, const char* fn);
void save_vtk(const PtsGraph& g, const char* fn);

//debugging flags for internal steps
void set_debug_mode(int i);
int get_debug_mode();

struct TicToc{
	//usage:
	//TicToc timer();
	//......
	//timer.fintoc();
	TicToc(bool start=true, const char* name = "Time duration");
	void tic();
	void toc();
	void report() const;
	void fintoc(){ toc(); report(); }
private:
	typedef std::chrono::high_resolution_clock TClock;
	typedef TClock::time_point TTimePoint;
	typedef std::chrono::duration<double> TDuration;
	const char* name;
	bool is_working;
	TTimePoint tp;
	TDuration dur;
};


#endif
