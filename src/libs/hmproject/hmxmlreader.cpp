#include "hmxmlreader.hpp"
#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#include <libxml/xpath.h>
#include <sstream>
#include <fstream>
#include <string.h>

using namespace HMXML;

std::vector<std::string> HMXML::split(const std::string& s, const std::string& delim){
	std::vector<std::string> ret;
	auto start = 0U;
	//auto end = s.find(delim);
	auto end = s.find_first_of(delim);
	while (end != std::string::npos){
		ret.push_back(s.substr(start, end - start));
		start = end + 1;
		end = s.find_first_of(delim, start);
	}
	ret.push_back(s.substr(start, end));
	
	std::vector<std::string> ret2;
	for (auto r: ret) if (r.size()>0) ret2.push_back(r);
	return ret2;
}

Reader::Reader(std::string fn):isroot(true){
	_doc=xmlParseFile(fn.c_str());
	if (_doc==NULL){
		throw std::runtime_error("Failed to parse xml document: "+fn);
	}
	_nd=xmlDocGetRootElement((xmlDoc*)_doc);
}
void Reader::Free(){
	if (isroot && _doc) {
		xmlFreeDoc((xmlDoc*)_doc);
	}
}
Reader Reader::create(std::string tag){
	Reader ret;
	xmlDoc* doc=xmlNewDoc((xmlChar*)"1.0");
	xmlNode* nd=xmlNewNode(NULL, (xmlChar*)tag.c_str());
	xmlDocSetRootElement(doc, nd);
	ret._doc = doc;
	ret._nd = nd;
	ret.isroot=true;
	return ret;
}

Reader Reader::copy_to_root(Reader& rd){
	Reader ret;
	ret._doc = xmlNewDoc((xmlChar*)"1.0");
	ret._nd = xmlCopyNode((xmlNode*)rd._nd, 1);
	xmlDocSetRootElement((xmlDoc*)ret._doc, (xmlNode*)ret._nd);
	ret.isroot=true;
	return ret;
}

std::string Reader::tostring(){
	xmlChar* s;
	int size;
	xmlDocDumpFormatMemoryEnc((xmlDoc*)_doc, &s, &size, "utf-8", 1);
	if (s==0) throw std::runtime_error("failed to write xml document");
	std::string ret((char*)s);
	xmlFree(s);
	size_t p1 = ret.find_first_not_of("\n\r\t\v ");
	size_t p2 = ret.find_last_not_of("\n\r\t\v ")+1;
	ret = ret.substr(p1, p2-p1);
	return ret;
}
void Reader::write(std::string filename){
	assert(isroot = true);
	xmlSaveFormatFileEnc(filename.c_str(), (xmlDoc*)_doc, "utf-8", 1);
}

Reader Reader::find_by_path(std::string path, bool required){
	xmlXPathContext *xpathCtx = xmlXPathNewContext((xmlDoc*)_doc);
	xpathCtx->node=(xmlNode*)_nd;
	xmlXPathObject *xpathObj = xmlXPathEvalExpression((xmlChar*)path.c_str(), xpathCtx);
	Reader ret;
	if (xpathObj!=NULL && xpathObj->nodesetval!=NULL && xpathObj->nodesetval->nodeNr!=0){
		for (int i=0; i<xpathObj->nodesetval->nodeNr; ++i){
			ret = Reader(*this, xpathObj->nodesetval->nodeTab[i]);
		}
	}
	xmlXPathFreeObject(xpathObj);
	xmlXPathFreeContext(xpathCtx);

	if (ret || !required) return ret;
	else throw XmlElementNotFound(path);
}

vector<Reader> Reader::findall_by_path(std::string path){
	xmlXPathContext *xpathCtx = xmlXPathNewContext((xmlDoc*)_doc);
	xpathCtx->node=(xmlNode*)_nd;
	xmlXPathObject *xpathObj = xmlXPathEvalExpression((xmlChar*)path.c_str(), xpathCtx);
	vector<Reader> ret;
	if (xpathObj!=NULL && xpathObj->nodesetval!=NULL && xpathObj->nodesetval->nodeNr!=0){
		for (int i=0; i<xpathObj->nodesetval->nodeNr; ++i){
			ret.push_back(Reader(*this, xpathObj->nodesetval->nodeTab[i]));
		}
	}
	xmlXPathFreeObject(xpathObj);
	xmlXPathFreeContext(xpathCtx);
	return ret;
}

bool Reader::value_string(std::string path, std::string& val, bool required){
	Reader fnd = find_by_path(path);
	if (!fnd){
		if (required) throw XmlElementNotFound(path);
		else return false;
	} else {
		val="";
		auto it=((xmlNode*)fnd._nd)->children;
		while (it!=NULL){
			if (it->type==3){
				auto cnt=xmlNodeListGetString((xmlDoc*)_doc, it, 0);
				val=(char*)cnt;
				xmlFree(cnt);
				break;
			} else it=it->next;
		}
		return true;
	}
}

std::string Reader::attribute(std::string path, std::string aname){
	Reader fnd = find_by_path(path);
	if (!fnd) return "";

	std::string ret;
	xmlAttr* attribute = ((xmlNode*)fnd._nd)->properties;
	while(attribute && attribute->name && attribute->children){
		xmlChar* value = xmlNodeListGetString((xmlDoc*)_doc, attribute->children, 1);
		std::string nm=(char*)attribute->name;
		if (nm == aname){
			ret=(char*)value;
			attribute = 0;
		} else attribute = attribute->next;
		xmlFree(value); 
	}
	return ret;
}
std::string Reader::tag_name(){
	return std::string((char*)((xmlNode*)_nd)->name);
}


bool Reader::value_float(std::string path, double& val, bool required){
	std::string s;
	if (!value_string(path, s)){
		if (required) throw XmlElementNotFound(path);
		else return false;
	} else {
		val = atof(s.c_str());
		return true;
	}
}

bool Reader::value_float_vec(std::string path, vector<double>& val, bool required){
	std::string s;
	if (!value_string(path, s)){
		if (required) throw XmlElementNotFound(path);
		else return false;
	}

	std::string attr = attribute(path, "tp");
	if (attr == "const" || attr == "constant"){
		double v = atof(s.c_str());
		for (auto& it: val) it=v;
	} else if (attr == "" || attr == "ascii"){
		vector<std::string> ss = split(s, " \n\r\t");
		if (val.size() < ss.size()) val.resize(ss.size());
		for (int i=0; i<ss.size(); ++i) val[i] = atof(ss[i].c_str());
	} else throw std::runtime_error("unknown vector data attribute");
	return true;
}

bool Reader::value_int(std::string path, int& val, bool required){
	std::string s;
	if (!value_string(path, s)){
		if (required) throw XmlElementNotFound(path);
		else return false;
	} else {
		val = atoi(s.c_str());
		return true;
	}
}

bool Reader::value_ulong(std::string path, unsigned long& val, bool required){
	std::string s;
	if (!value_string(path, s)){
		if (required) throw XmlElementNotFound(path);
		else return false;
	} else {
		val = std::stoul(s.c_str());
		return true;
	}
}

void Reader::new_attribute(std::string key, std::string a){
	xmlSetProp((xmlNode*)_nd, (xmlChar*)key.c_str(), (xmlChar*)a.c_str());
}
void Reader::delete_content(){
	auto it=((xmlNode*)_nd)->children;
	while (it!=NULL){
		if (it->type==3){
			auto it2=it->next;
			xmlUnlinkNode(it);
			xmlFreeNode(it);
			it=it2;
		} else it=it->next;
	}
}
	
void Reader::unlink_node(){
	if (_nd != nullptr){
		xmlUnlinkNode((xmlNode*)_nd);
		xmlFreeNode((xmlNode*)_nd);
		_nd = nullptr; _doc = nullptr;
	}
}

void Reader::set_content(const std::string& data){
	delete_content();
	auto tnd=xmlNewText((xmlChar*)data.c_str());
	xmlAddChild((xmlNode*)_nd, tnd);
}

Reader Reader::new_child(std::string tag){
	xmlNode* nd=(xmlNode*)_nd;
	xmlNode* cnode=xmlNewChild(nd, NULL, (xmlChar*)tag.c_str(), NULL);
	return Reader(*this, cnode);
}
Reader Reader::provide_child(std::string tag){
	auto fnd = find_by_path(tag);
	if (!fnd) fnd = new_child(tag);
	return fnd;
}

// ======================== AReader
namespace{
template<class A>
shared_ptr<std::ostringstream> getstream(){
	return shared_ptr<std::ostringstream>(new std::ostringstream());
}
template<> shared_ptr<std::ostringstream> getstream<double>(){
	shared_ptr<std::ostringstream> ret(new std::ostringstream());
	ret->precision(16);
	return ret;
}
template<> shared_ptr<std::ostringstream> getstream<float>(){
	shared_ptr<std::ostringstream> ret(new std::ostringstream());
	ret->precision(8);
	return ret;
}
template<class A>
std::string tostr(const vector<A>& v){
	auto s=getstream<A>();
	for (auto x: v) (*s)<<x<<' ';
	return s->str();
};
template<> std::string tostr<char>(const vector<char>& v){
	auto s=getstream<char>();
	for (auto x: v) (*s)<<(int)x<<' ';
	return s->str();
};
template<class A>
std::string vtostr(const vector<vector<A>>& v, bool isvariable){
	if (!isvariable){
		std::vector<A> alldata;
		for (auto& it: v) std::copy(it.begin(), it.end(), std::back_inserter(alldata));
		return tostr(alldata);
	} else {
		std::ostringstream s;
		for (auto& it: v){
			s<<(int)it.size()<<" ";
			s<<tostr(it);
		}
		return s.str();
	}
}

template<class A>
std::string get_type_id(){
	throw std::runtime_error("Unsupported type");
};
template<> std::string get_type_id<int>(){return "int";}
template<> std::string get_type_id<char>(){return "char";}
template<> std::string get_type_id<double>(){return "double";}
template<> std::string get_type_id<float>(){return "float";}

template<class A>
void add_field(HMXML::Reader& reader, const std::vector<A>& data, bool binary, vector<char>& buffer){
	if (data.size() == 0) return;
	reader.new_attribute("type", get_type_id<A>());
	if (!binary){
		reader.new_attribute("format", "ascii");
		reader.set_content(tostr(data));
	} else {
		reader.new_attribute("format", "binary");
		size_t sz= data.size()*sizeof(A);
		reader.new_child("START").set_content(std::to_string(buffer.size()));
		buffer.resize(buffer.size()+sz);
		std::copy((char*)(&data[0]), (char*)(&data[0])+sz, buffer.end()-sz);
	}
}
template<class A>
void add_vfield(HMXML::Reader& reader, const vector<vector<A>>& data, bool binary, vector<char>& buffer){
	if (data.size() == 0) return;
	reader.new_attribute("type", get_type_id<A>());
	int vecsz=data[0].size();
	for (int i=1; i<data.size(); ++i) if (data[i].size() != vecsz){
		vecsz = -1; break;
	}
	int totalsz=0;
	for (int i=0; i<data.size(); ++i) totalsz+=data[i].size();

	if (vecsz == -1) reader.new_attribute("dim", "variable");
	else if (vecsz != 1) reader.new_attribute("dim", std::to_string(vecsz));
	if (!binary){
		reader.new_attribute("format", "ascii");
		reader.set_content(vtostr(data, vecsz==-1));
	} else {
		reader.new_attribute("format", "binary");
		size_t sz= totalsz*sizeof(A);
		if (vecsz == -1) sz += data.size()*sizeof(unsigned int);
		reader.new_child("START").set_content(std::to_string(buffer.size()));
		buffer.resize(buffer.size()+sz);
		auto bufiter = buffer.end()-sz;
		if (vecsz != -1){
			size_t itsz = vecsz*sizeof(A);
			for (auto& it: data){
				std::copy((char*)(&it[0]), (char*)(&it[0])+itsz, bufiter);
				bufiter+=itsz;
			}
		} else {
			for (auto& it: data){
				unsigned int itsz = it.size();
				std::copy((char*)(&itsz), (char*)(&itsz)+sizeof(unsigned int), bufiter);
				bufiter+=sizeof(unsigned int);
				itsz*=sizeof(A);
				std::copy((char*)(&it[0]), (char*)(&it[0])+itsz, bufiter);
				bufiter+=itsz;
			}
		}
	}
}
}//namespace

void ReaderA::set_num_content(const std::vector<char>& data, Reader& subnode, bool binary){
	add_field(subnode, data, binary, buffer);
}
void ReaderA::set_num_content(const std::vector<int>& data, Reader& subnode, bool binary){
	add_field(subnode, data, binary, buffer);
}
void ReaderA::set_num_content(const std::vector<float>& data, Reader& subnode, bool binary){
	add_field(subnode, data, binary, buffer);
}
void ReaderA::set_num_content(const std::vector<double>& data, Reader& subnode, bool binary){
	add_field(subnode, data, binary, buffer);
}

void ReaderA::set_num_content(const std::vector<std::vector<char>>& data, Reader& subnode, bool binary){
	add_vfield(subnode, data, binary, buffer);
}
void ReaderA::set_num_content(const std::vector<std::vector<int>>& data, Reader& subnode, bool binary){
	add_vfield(subnode, data, binary, buffer);
}
void ReaderA::set_num_content(const std::vector<std::vector<float>>& data, Reader& subnode, bool binary){
	add_vfield(subnode, data, binary, buffer);
}
void ReaderA::set_num_content(const std::vector<std::vector<double>>& data, Reader& subnode, bool binary){
	add_vfield(subnode, data, binary, buffer);
}

ReaderA::ReaderA(Reader& r):Reader(r){
	assert(r.isroot);
}

ReaderA ReaderA::create(std::string tag){
	Reader r = Reader::create(tag);
	return ReaderA(r);
}
ReaderA* ReaderA::pcreate(std::string tag){
	Reader r = Reader::create(tag);
	return new ReaderA(r);
}

void ReaderA::write(std::string filename){
	std::ofstream ofile(filename, std::ios::out | std::ios::binary);
	std::string xmlstring = tostring();
	while (xmlstring.size()!=0 && xmlstring.back()!='>') xmlstring.pop_back();
	assert(xmlstring.size() > 0);
	ofile<<xmlstring;
	if (buffer.size()>0){
		ofile.write(&buffer[0], buffer.size());
	}
}

ReaderA::ReaderA(std::string fn, std::string ending_tag){
	std::ifstream fs(fn);
	if (!fs.is_open()) throw std::runtime_error("file "+fn+" was not found");
	auto buf = fs.rdbuf();
	size_t sz = buf->pubseekoff(0, fs.end, fs.in);
	buf->pubseekpos(0, fs.in);
	vector<char> bufc(sz);
	buf->sgetn(&bufc[0], sz);
	auto fnd = strstr(&bufc[0], ending_tag.c_str());
	if (fnd == NULL) throw std::runtime_error("ending tag "+ending_tag+" was not found");
	fnd += ending_tag.size();

	size_t fndpos = (fnd-&bufc[0]);
	size_t endpos = bufc.size();

	_doc = xmlParseMemory(&bufc[0], fndpos);
	if (_doc == NULL) std::runtime_error("Error reading xml");
	_nd=xmlDocGetRootElement((xmlDoc*)_doc);
	isroot=true;

	size_t bufsize = bufc.size()-fndpos;
	buffer.resize(bufsize);
	std::copy(bufc.begin()+fndpos, bufc.end(), buffer.begin());

	if (buffer.size() == 1 && buffer[0]=='\n') buffer.clear();
	if (buffer.size() == 2 && buffer[0]=='\n' && buffer[1]=='\r') buffer.clear();
}

ReaderA::TNumContent::TNumContent(Reader& subnode, int num, const vector<char>& buffer){
	if (num==0) return;

	//read attributes
	std::string type_str = subnode.attribute(".", "type");
	std::string dim_str = subnode.attribute(".", "dim");
	std::string format_str=subnode.attribute(".", "format");

	//process attributes
	if (dim_str=="") dim_str="1";
	if (dim_str=="variable") dim_str="-1";
	int dim = atoi(dim_str.c_str());

	//read data
	if (format_str == "ascii"){
		std::string content;
		subnode.value_string(".", content, true);
		fill_data_ascii(content, dim, type_str, num);
	} else if (format_str == "binary"){
		unsigned long pos;
		subnode.value_ulong("START", pos, true);
		fill_data_bin(&buffer[0]+pos, buffer.size()-pos, dim, type_str, num);
	} else throw std::runtime_error("unknown xml field format "+format_str);
}

ReaderA::TNumContent ReaderA::read_num_content(Reader& subnode, int num){
	return ReaderA::TNumContent(subnode, num, buffer);
}

namespace{
template<class A> void toA(const std::string& s, A& val){}; 
template<> void toA<float>(const std::string& s, float& val){val=(float)atof(s.c_str());}
template<> void toA<double>(const std::string& s, double& val){val=atof(s.c_str());}
template<> void toA<int>(const std::string& s, int& val){val=atoi(s.c_str());}
template<> void toA<char>(const std::string& s, char& val){val=(char)atoi(s.c_str());}

template<class A>
void fill_vector(vector<A>& ret, vector<std::string>::const_iterator start, vector<std::string>::const_iterator end){
	ret.resize(end-start);
	auto rit = ret.begin();
	for (auto it=start; it!=end; ++it) toA(*it, *rit++);
}
template<class A>
void fill_vector_bin(const char* start, size_t size, vector<A>& vec, int vecn){
	vec.resize(vecn);
	size_t asize = sizeof(A)*vecn;
	if (asize > size) throw std::runtime_error("file binary buffer overflow");
	memcpy((void*)(&vec[0]), start, asize);
}

template<class A>
void split_vector(vector<A>& input, vector<vector<A>>& output, int dim){
	output = vector<vector<A>>(input.size()/dim, vector<A>(dim));
	auto it=input.begin();
	for (int i=0; i<output.size(); ++i){
		for (int j=0; j<dim; ++j){
			output[i][j]=*it++;
		}
	}
	input.clear();
}
template<class A>
void fill_vvector(const vector<std::string>& s, vector<vector<A>>& ret){
	auto it=s.begin();
	while (it!=s.end()){
		int dim = atoi((*it++).c_str());
		if (it-s.begin()+dim>s.size())
			throw std::runtime_error("corrupted variable vector data");
		ret.emplace_back();
		fill_vector<A>(ret.back(), it, it+dim);
		it+=dim;
	}
}
template<class A>
void fill_vvector_bin(const char* start, size_t size, vector<vector<A>>& ret, int vecn){
	ret.resize(vecn);
	const char* iter=start;
	long leftsize=size;
	auto decrease_leftsize = [&leftsize](size_t num){
		leftsize-=num;
		if (leftsize<0) throw std::runtime_error("file binary buffer overflow");
	};
	for (int i=0; i<vecn; ++i){
		//dim
		unsigned int dim;
		decrease_leftsize(sizeof(unsigned int));
		memcpy(&dim, iter, sizeof(unsigned int));
		iter+=sizeof(unsigned int);
		//data
		fill_vector_bin(iter, leftsize, ret[i], dim);
		decrease_leftsize(dim*sizeof(A));
		iter+=dim*sizeof(A);
	}
}

}
void ReaderA::TNumContent::fill_data_ascii(const std::string& data, int dim, std::string typestr, int num){
	int tp = code(typestr);
	if (dim == 1){
		std::vector<std::string> ss = split(data, "\n\r\v\t ");
		switch (tp){
			case RINT: fill_vector(intvec, ss.begin(), ss.end()); break;
			case RCHR: fill_vector(chrvec, ss.begin(), ss.end()); break;
			case RDBL: fill_vector(dblvec, ss.begin(), ss.end()); break;
			case RFLT: fill_vector(fltvec, ss.begin(), ss.end()); break;
		};
	} else if (dim!=-1){
		fill_data_ascii(data, 1, typestr, num*dim);
		switch (tp){
			case RINT: split_vector(intvec, intvecvec, dim); break;
			case RCHR: split_vector(chrvec, chrvecvec, dim); break;
			case RDBL: split_vector(dblvec, dblvecvec, dim); break;
			case RFLT: split_vector(fltvec, fltvecvec, dim); break;
		}
	} else{
		std::vector<std::string> ss = split(data, "\n\r\v\t ");
		switch (tp){
			case RINT: fill_vvector(ss, intvecvec); break;
			case RCHR: fill_vvector(ss, chrvecvec); break;
			case RDBL: fill_vvector(ss, dblvecvec); break;
			case RFLT: fill_vvector(ss, fltvecvec); break;
		}
	}
}
void ReaderA::TNumContent::fill_data_bin(const char* start, size_t size, int dim, std::string typestr, int num){
	int tp = code(typestr);
	if (dim == 1){
		switch (tp){
			case RINT: fill_vector_bin(start, size, intvec, num); break;
			case RFLT: fill_vector_bin(start, size, fltvec, num); break;
			case RDBL: fill_vector_bin(start, size, dblvec, num); break;
			case RCHR: fill_vector_bin(start, size, chrvec, num); break;
		}
	} else if (dim!=-1){
		fill_data_bin(start, size, 1, typestr, num*dim);
		switch (tp){
			case RINT: split_vector(intvec, intvecvec, dim); break;
			case RCHR: split_vector(chrvec, chrvecvec, dim); break;
			case RDBL: split_vector(dblvec, dblvecvec, dim); break;
			case RFLT: split_vector(fltvec, fltvecvec, dim); break;
		}
	} else {
		switch (tp){
			case RINT: fill_vvector_bin(start, size, intvecvec, num); break;
			case RCHR: fill_vvector_bin(start, size, chrvecvec, num); break;
			case RDBL: fill_vvector_bin(start, size, dblvecvec, num); break;
			case RFLT: fill_vvector_bin(start, size, fltvecvec, num); break;
		}
	}
	if (N()!=num) throw std::runtime_error("invalid ascii numerical field size given");
}
int ReaderA::TNumContent::N() const{
	if (intvec.size()>0) return intvec.size();
	if (fltvec.size()>0) return fltvec.size();
	if (dblvec.size()>0) return dblvec.size();
	if (chrvec.size()>0) return chrvec.size();
	if (intvecvec.size()>0) return intvecvec.size();
	if (fltvecvec.size()>0) return fltvecvec.size();
	if (dblvecvec.size()>0) return dblvecvec.size();
	if (chrvecvec.size()>0) return chrvecvec.size();
	return 0;
}
int ReaderA::TNumContent::code(std::string s) const{
	if (s=="int") return RINT;
	if (s=="float") return RFLT;
	if (s=="double") return RDBL;
	if (s=="char") return RCHR;
	throw std::runtime_error("unknown data type "+s);
}


namespace{
//This object is used to invoke xmlCleanupParser() at the end of the library utilization
//in order to prevent memory leaks.
struct CleanUp{
	~CleanUp(){ xmlCleanupParser(); }
} cleanup;
}//namespace
