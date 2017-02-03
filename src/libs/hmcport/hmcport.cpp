#include "hmcport.h"
#include "hmtesting.hpp"
#include <iostream>

int free_int_array(int* a){
	try{
		if (a != NULL) delete[] a;
		return HMSUCCESS;
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}
int free_double_array(double* a){
	try{
		if (a != NULL) delete[] a;
		return HMSUCCESS;
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}
int free_char_array(char* a){
	try{
		if (a != NULL) delete[] a;
		return HMSUCCESS;
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}
int free_voidp_array(void** a){
	try{
		if (a != NULL) delete[] a;
		return HMSUCCESS;
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}

int get_ascii_file_hash(const char* fn, size_t* ret){
	try{
		*ret = HMTesting::calculate_file_hash(fn);
		return HMSUCCESS;
	} catch (std::exception& e){
		std::cout<<e.what()<<std::endl;
		return HMERROR;
	}
}

const char* get_boundary_name(const BoundaryNamesStruct* bnames, int index){
	for (int i=0; i<bnames->n; ++i){
		if (bnames->index[i] == index){
			return bnames->name[i];
		}
	}
	std::string ret = "boundary"+std::to_string(index);
	return ret.c_str();
}
