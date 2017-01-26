#ifndef FLUENT_EXPORT_GRID2D_HPP
#define FLUENT_EXPORT_GRID2D_HPP

#include "primitives2d.hpp"

namespace HM2D{ namespace Export{ 
typedef std::function<std::string(int)> BNamesFun;

struct PeriodicData{
	void add_data(int bt1, int bt2, bool is_reversed=true){
		b1.push_back(bt1); b2.push_back(bt2); isrev.push_back(is_reversed); }
	void clear(){ b1.clear(); b2.clear(); isrev.clear(); }
	std::vector<int> assemble(GridData& g, const vector<int>&) const;
	int size() const { return b1.size(); }
	std::vector<int> b1;
	std::vector<int> b2;
	std::vector<bool> isrev;
};

//bndindex.size() == total number of grid edges
void GridMSH(const GridData& g, std::string fn, vector<int> bndindex);

void GridMSH(const GridData& g, std::string fn, vector<int> bndindex, BNamesFun bnames);

void GridMSH(const GridData& g, std::string fn, vector<int> bndindex, PeriodicData pd);

void GridMSH(const GridData& g, std::string fn, vector<int> bndindex, BNamesFun bnames, PeriodicData pd);

}}

#endif
