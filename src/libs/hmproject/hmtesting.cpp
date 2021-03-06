#include "hmtesting.hpp"
#include <fstream>
#include <sstream>
#include <iostream>

namespace{
int FAILED_CHECKS = 0;

std::vector<std::string> split(const std::string& s, const std::string& delim){
	std::vector<std::string> ret;
	auto start = 0U;
	auto end = s.find(delim);
	while (end != std::string::npos){
		ret.push_back(s.substr(start, end - start));
		start = end + delim.length();
		end = s.find(delim, start);
		//std::cout << s.substr(start, end);
	}
	ret.push_back(s.substr(start, end));
	return ret;
}

std::string remove_all_chars(const std::string& s, char c){
	std::vector<std::string> s2 = split(s, std::string(1, c));
	std::string ret;
	for (auto& s: s2){
		ret += s;
	}
	return ret;
}

std::string read_file_to_string(std::string fn){
	std::ifstream t(fn);
	std::stringstream buffer;
	buffer << t.rdbuf();
	return buffer.str();
}

}

std::string HMTesting::convert_floats(const std::string& inp, int rnd){
	std::string ret;
	std::string nums = "0123456789-+eE";

	double rmult = 1;
	for (int i=0; i<rnd; ++i) rmult*=10;
	auto rndfun = [&rmult](double& x){
		x = round(x*rmult)/rmult;
		if (ISZERO(x)) x = 0;  //to avoid '-0'
	};
	
	vector<std::string> sp = split(inp, ".");

	vector<std::string> mrg;
	for (int i=0; i<sp.size()-1; ++i){
		auto& s1 = sp[i];
		auto& s2 = sp[i+1];
		size_t p1 = s1.find_last_not_of(nums);
		if (p1 == std::string::npos) p1 = 0;
		else (p1 += 1);
		size_t p2 = s2.find_first_not_of(nums);
		if (p2 == std::string::npos) p2 = s2.size();

		std::string snum = s1.substr(p1, s1.size()) + "." + s2.substr(0, p2);
		double a;
		try{
			a = std::stod(snum);
			rndfun(a);
		} catch (...){
			mrg.push_back(".");
			std::cout<<"error at "<<snum<<std::endl;
			continue;
		}
		s1.resize(p1);
		s2 = s2.substr(p2, s2.size());

		std::stringstream sstream;
		sstream<<a;
		mrg.push_back(sstream.str());
	}
	for (int i=0; i<sp.size()-1; ++i){
		ret += sp[i];
		ret += mrg[i];
	}
	ret += sp.back();
	return ret;
}

int HMTesting::failed_num(){
	return FAILED_CHECKS;
}
void HMTesting::check_final_report(){
	if (FAILED_CHECKS ==1){
		std::cout<<FAILED_CHECKS<<" test failed <<<<<<<<<<<<<<<<<<<"<<std::endl;
	} else if (FAILED_CHECKS > 1) {
		std::cout<<FAILED_CHECKS<<" tests failed <<<<<<<<<<<<<<<<<<<"<<std::endl;
	} else {
		std::cout<<"All tests passed"<<std::endl;
	}
}

void HMTesting::add_check(bool ex, std::string info){
	if (info.size()==0){
		std::cout<<"\tunknown check: ";
	} else{
		std::cout<<"\t"<<info;
	}
	if (ex){
		std::cout<<": True"<<std::endl;
	} else {
		++FAILED_CHECKS;
		std::cout<<": False <<<<<<<<<<<<<<<<<<<"<<std::endl;
	}
};

size_t HMTesting::calculate_file_hash(std::string fn){
	auto str = read_file_to_string(fn);
	str = convert_floats(str, 2);
	str = remove_all_chars(str, '\r');
	return std::hash<std::string>()(str);
}

void HMTesting::add_file_check(size_t hash, std::string fn, std::string info){
	size_t h1 = calculate_file_hash(fn);
	int fc = FAILED_CHECKS;
	add_check(h1 == hash, info);
	if (fc != FAILED_CHECKS){
		std::cout<<"\t\tcurrent file hash = "<<h1<<'U'<<std::endl;
		std::cout<<"\t\tinput file hash   = "<<hash<<'U'<<std::endl;
	}
}

void HMTesting::add_file_check(std::string f1, std::string f2, std::string info){
	size_t h1 = calculate_file_hash(f1);
	size_t h2 = calculate_file_hash(f2);
	add_check(h1 == h2, info);
}
size_t HMTesting::calculate_file_size(std::string fn){
	std::ifstream in(fn, std::ifstream::ate | std::ifstream::binary);
	return in.tellg(); 
}

void HMTesting::add_filesize_check(size_t sz, std::string fn, std::string info){
	size_t h1 = calculate_file_size(fn);
	int fc = FAILED_CHECKS;
	add_check(h1 == sz, info);
	if (fc != FAILED_CHECKS){
		std::cout<<"\t\tcurrent file size = "<<h1<<std::endl;
		std::cout<<"\t\tinput file size   = "<<sz<<std::endl;
	}
}
