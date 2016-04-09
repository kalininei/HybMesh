#ifndef HYBMESH_TESTING_HPP
#define HYBMESH_TESTING_HPP
#include "hmproject.h"

namespace HMTesting{
	//reads all float numbers in input string, rounds them
	//to 'round' digits after a point and writes them back
	std::string convert_floats(const std::string& inp, int round);

	//check condition
	void add_check(bool ex, std::string info);
	
	//compares file with given hash value with floats rounding
	void add_file_check(size_t hash, std::string fn, std::string info);

	//calculate ascii file content hash with float rounding
	size_t calculate_file_hash(std::string fn);

	//returns number of failed checks;
	int failed_num();

	//prints number of failed tests or ok
	void check_final_report();
}

#endif
