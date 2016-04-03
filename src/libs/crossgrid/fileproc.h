#ifndef CROSSGRID_FILEPROC_H
#define CROSSGRID_FILEPROC_H
#include <chrono>
#include "bgeom2d.h"
#include "grid.h"
#include "wireframegrid.h"

// ---------------- debugging procedures
void save_vtk(const GridGeom* g, const char* fn);
void save_vtk(const GridGeom* g, const vector<double>& pdata, const char* fn);
void save_vtk(const PContour* c, const char* fn);
void save_vtk(const vector<PContour>& c, const char* fn);
void save_vtk(const vector<PContour>& c, const vector<double>& data, const char* fn);
void save_vtk(const PtsGraph* g, const char* fn);

void save_vtk(const ContoursCollection& c, const char* fn);
void save_vtk(const ContoursCollection* c, const char* fn);
void save_vtk(const ContoursCollection& c, const vector<double>& pdata, const char* fn);
void save_vtk(const GridGeom& g, const char* fn);
void save_vtk(const GridGeom& g, const vector<double>& pdata, const char* fn);
void save_vtk(const PtsGraph& g, const char* fn);

#endif
