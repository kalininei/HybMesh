#ifndef HYBMESH_XML_READER_HPP
#define HYBMESH_XML_READER_HPP
#include "hmproject.h"

namespace HMXML{

std::vector<std::string> split(const std::string& s, const std::string& delim);

struct Reader{
	Reader(std::string fn);
	Reader(Reader& parent, void* node):isroot(false), _doc(parent._doc), _nd(node){};
	Reader():isroot(false), _doc(nullptr), _nd(nullptr){};
	virtual ~Reader(){}
	//we cannot rely on standrard destructor because in that case we'll have
	//problems with copying of root elements. So user should use explicit free for a root
	//node document.
	void Free();
	std::string tostring();
	virtual void write(std::string filename);

	//===== Read procedures
	//returns true if data was read
	bool value_string(std::string path, std::string& val, bool required=false);
	bool value_float(std::string path, double& val, bool required=false);
	bool value_float_vec(std::string path, vector<double>& val, bool required=false);
	bool value_int(std::string path, int& val, bool required=false);
	bool value_ulong(std::string path, unsigned long& val, bool required=false);
	std::string attribute(std::string path, std::string aname);
	std::string tag_name();

	//===== Modify procedures
	Reader new_child(std::string tag);
	Reader provide_child(std::string tag);
	void new_attribute(std::string key, std::string);
	void set_content(const std::string& data);
	void delete_content();
	void unlink_node();

	//===== xpath find procedures
	Reader find_by_path(std::string path, bool required=false);
	vector<Reader> findall_by_path(std::string path);

	operator bool() const {return _nd != nullptr;}
	static Reader create(std::string tag);

	bool isroot;
	void *_doc, *_nd;
};

//reader with additional binary buffer, which will be written to end of output file
//this class should represent only the root node of xml document
struct ReaderA: public Reader{
	std::vector<char> buffer;
	ReaderA():Reader(){isroot=true;}
	ReaderA(Reader& r);
	ReaderA(std::string filename, std::string ending_tag);

	static ReaderA create(std::string tag);
	static ReaderA* pcreate(std::string tag);

	void set_num_content(const std::vector<char>& data, Reader& subnode, bool binary);
	void set_num_content(const std::vector<int>& data, Reader& subnode, bool binary);
	void set_num_content(const std::vector<float>& data, Reader& subnode, bool binary);
	void set_num_content(const std::vector<double>& data, Reader& subnode, bool binary);

	void set_num_content(const std::vector<std::vector<char>>& data, Reader& subnode, bool binary);
	void set_num_content(const std::vector<std::vector<int>>& data, Reader& subnode, bool binary);
	void set_num_content(const std::vector<std::vector<float>>& data, Reader& subnode, bool binary);
	void set_num_content(const std::vector<std::vector<double>>& data, Reader& subnode, bool binary);

	struct TNumContent;
	TNumContent read_num_content(Reader& subnode, int num);

	void write(std::string filename) override;
};

struct ReaderA::TNumContent{
	static const int RCHR=1;
	static const int RINT=2;
	static const int RFLT=3;
	static const int RDBL=4;
	TNumContent(Reader& subnode, int num, const vector<char>& buffer);

	//only one of those fields are filled after read_num_content execution
	std::vector<char> chrvec;
	std::vector<int> intvec;
	std::vector<float> fltvec;
	std::vector<double> dblvec;
	std::vector<std::vector<char>> chrvecvec;
	std::vector<std::vector<int>> intvecvec;
	std::vector<std::vector<float>> fltvecvec;
	std::vector<std::vector<double>> dblvecvec;
	int N() const;
	int code(std::string s) const;


	//return vector by types: int, vector<int> etc
	template<class A>
	std::vector<A>& vec(){
		force_v<A>();
		return _vec<A>();
	}
	//forces vector<vector<tp>> return even if dim=1
	template<class A>
	std::vector<std::vector<A>>& vecvec(){
		force_vv<A>(); 
		return _vecvec<A>();
	}
	template<class From, class To>    //From is one of int, float, char, double
	std::vector<To> convert_data();   //To is a simple class
	template<class From, class To>
	std::vector<vector<To>> convert_vdata();  //To is vector<...> class

private:
	void fill_data_ascii(const std::string& data, int dim, std::string typestr, int num);
	void fill_data_bin(const char* start, size_t size, int dim, std::string type_str, int num);
	template<class A> void force_vv();
	template<class A> void force_v();
	template<class A> vector<A>& _vec();
	template<class A> vector<vector<A>>& _vecvec();
	template<class A>
	static void rewrite_vecs_tovv(const vector<A>& v1, vector<vector<A>>& v2){
		if (v1.size()>0 && v2.size() == 0){
			v2.resize(v1.size(), vector<A>(1));
			for (int i=0; i<v2.size(); ++i) v2[i][0] = v1[i];
		}
	};
	template<class A>
	static void rewrite_vecs_tov(const vector<vector<A>>& v1, vector<A>& v2){
		if (v1.size()>0 && v2.size() == 0){
			int totdim = 0;
			for (auto& it: v1) totdim+= it.size();
			v2.resize(totdim);
			auto it=v2.begin();
			for (int i=0; i<v1.size(); ++i){
				std::copy(v1[i].begin(), v1[i].end(), it);
				it+=v1[i].size();
			}
		}
	};
};

template<> inline std::vector<char>& ReaderA::TNumContent::_vec<char>(){ return chrvec; }
template<> inline std::vector<int>& ReaderA::TNumContent::_vec(){ return intvec; }
template<> inline std::vector<float>& ReaderA::TNumContent::_vec(){ return fltvec; }
template<> inline std::vector<double>& ReaderA::TNumContent::_vec(){ return dblvec; }

template<> inline std::vector<vector<char>>& ReaderA::TNumContent::_vecvec(){ return chrvecvec; }
template<> inline std::vector<vector<int>>& ReaderA::TNumContent::_vecvec(){ return intvecvec; }
template<> inline std::vector<vector<float>>& ReaderA::TNumContent::_vecvec(){ return fltvecvec; }
template<> inline std::vector<vector<double>>& ReaderA::TNumContent::_vecvec(){ return dblvecvec; }

template<> inline void ReaderA::TNumContent::force_vv<int>(){ rewrite_vecs_tovv(intvec, intvecvec); }
template<> inline void ReaderA::TNumContent::force_vv<char>(){ rewrite_vecs_tovv(chrvec, chrvecvec); }
template<> inline void ReaderA::TNumContent::force_vv<double>(){ rewrite_vecs_tovv(dblvec, dblvecvec); }
template<> inline void ReaderA::TNumContent::force_vv<float>(){ rewrite_vecs_tovv(fltvec, fltvecvec); }

template<> inline void ReaderA::TNumContent::force_v<int>(){ rewrite_vecs_tov(intvecvec, intvec); }
template<> inline void ReaderA::TNumContent::force_v<char>(){ rewrite_vecs_tov(chrvecvec, chrvec); }
template<> inline void ReaderA::TNumContent::force_v<double>(){ rewrite_vecs_tov(dblvecvec, dblvec); }
template<> inline void ReaderA::TNumContent::force_v<float>(){ rewrite_vecs_tov(fltvecvec, fltvec); }


template<class From, class To>
std::vector<To> ReaderA::TNumContent::convert_data(){
	auto& vv = vec<From>();
	std::vector<To> ret(vv.size());
	std::copy(vv.begin(), vv.end(), ret.begin());
	return ret;
}
template<class From, class To>
std::vector<vector<To>> ReaderA::TNumContent::convert_vdata(){
	auto& vv = vecvec<From>();
	std::vector<vector<To>> ret(vv.size());
	for (int i=0; i<ret.size(); ++i){
		std::copy(vv[i].begin(), vv[i].end(), std::back_inserter(ret[i]));
	}
	return ret;
}




};



#endif
