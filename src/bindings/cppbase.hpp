#ifndef HYBMESH_CPP_FRONTEND_HPP
#define HYBMESH_CPP_FRONTEND_HPP

#include <vector>
#include <string>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <cstdlib>
#include <functional>

#define DEFAULT_HYBMESH_EXE_PATH "hybmesh"  // >>$EXEPATH

class Hybmesh{
public:
	struct Point2;
	struct Point3;
	struct Contour2D;
	struct Grid2D;
	struct Surface3D;
	struct Grid3D;
	typedef std::vector<char> VecByte;
	struct ERuntimeError: public std::runtime_error{
		ERuntimeError(const char* msg) noexcept: std::runtime_error(
			std::string("Hybmesh runtime error: ") +
			std::string(msg)){};
	};
	struct EUserInterrupt: public std::runtime_error{
		EUserInterrupt() noexcept: std::runtime_error(
			std::string("Interrupted by user")){};
	};
	struct Point2{
		double x, y;
		Point2(double x=0, double y=0): x(x), y(y){}
		static Point2 None(){
			return Point2(std::numeric_limits<double>::max(),
					std::numeric_limits<double>::max());
		}
		bool isnone() const{
			return x == std::numeric_limits<double>::max() &&
			       y == std::numeric_limits<double>::max();
		}
	};
	struct Point3{
		double x, y, z;
		Point3(double x=0, double y=0, double z=0): x(x), y(y), z(z){}
		static Point3 None(){
			return Point3(std::numeric_limits<double>::max(),
					std::numeric_limits<double>::max(),
					std::numeric_limits<double>::max());
		}
		bool isnone() const{
			return x == std::numeric_limits<double>::max() &&
			       y == std::numeric_limits<double>::max() &&
			       z == std::numeric_limits<double>::max();
		}
	};
private:
	class Worker{
		//data
		int sig_read, sig_write, data_read, data_write;
		intptr_t childid;
		//low level platform specific procedures
		void require_connection(const char* path);
		void get_signal(char* sig);
		void get_data(int* sz, char** data);
		void send_signal(char sig);
		void send_data(int sz, const char* data);
		void break_connection();
		//server communication
		void _send_command(const std::string& func, const std::string& com);
		char _wait_for_signal();
		VecByte _read_buffer();
	public:
		//callback
		std::function<int(const std::string&, const std::string&,
		                  double, double)> callback;
		//constructor
		Worker(const char* hybmeshpath){
			if (hybmeshpath == NULL) require_connection(
					DEFAULT_HYBMESH_EXE_PATH);
			else require_connection(hybmeshpath);
			callback = [](const std::string&, const std::string&,
			              double, double)->int{ return 0; };
		};
		~Worker(){
			break_connection();
		};

		//command interface
		VecByte _apply_command(const std::string& func, const std::string& com);
		void _apply_callback(const VecByte& buf);
		// string -> cpp type
		int _to_int(const std::string&);
		double _to_double(const std::string&);
		Point2 _to_point(const std::string&);
		Grid2D _to_grid(const std::string&);
		Contour2D _to_cont(const std::string&);
		Surface3D _to_surface(const std::string&);
		Grid3D _to_grid3(const std::string&);
		std::vector<std::string> __parse_vecstring(const std::string&);
		std::vector<int> _to_vecint(const std::string&);
		std::vector<Contour2D> _to_veccont(const std::string&);
		std::vector<Grid2D> _to_vecgrid(const std::string&);
		std::vector<Surface3D> _to_vecsurface(const std::string&);
		std::vector<Grid3D> _to_vecgrid3(const std::string&);
		template<class Obj> Obj _to_object(const std::string&);
		template<class Obj> std::vector<Obj> _to_vecobject(const std::string&);
		// cpp type -> string
		std::string _tos_bool(bool);
		std::string _tos_int(int);
		std::string _tos_double(double);
		std::string _tos_point(Point2);
		std::string _tos_point3(Point3);
		std::string _tos_vecbyte(const VecByte& comstr);
		std::string _tos_vecint(const std::vector<int>&);
		std::string _tos_vecdouble(const std::vector<double>&);
		std::string _tos_vecstring(const std::vector<std::string>&);
		std::string _tos_vecpoint(const std::vector<Point2>&);
		std::string _tos_string(const std::string&);
		template<class Obj> std::string _tos_object(const Obj&);
		template<class Obj> std::string _tos_vecobject(const std::vector<Obj>&);
		// vecbyte -> cpp type
		std::vector<std::pair<int, double> > _to_vec_int_double_raw(const VecByte&);
		std::vector<double> _to_vecdouble_raw(const VecByte&);
		std::vector<int> _to_vecint_raw(const VecByte&);
	};

	Worker _worker;
	Worker* worker;

	// copy is forbidden
	Hybmesh(const Hybmesh&);
public:
	// =============== construction/destruction
	Hybmesh(const char* hybmeshpath=NULL): _worker(hybmeshpath) {
		worker = &_worker;
	}
	// ================ callback
	void assign_callback(std::function<int(const std::string&, const std::string&,
		                               double, double)> func){
		worker->callback = func;
	}
	void reset_callback(){
		worker->callback = [](const std::string&, const std::string&,
		                      double, double)->int { return 0; };
	}
	// ================ Geometrical objects
	struct Object{
		std::string sid;
		Worker* worker;
		Object(std::string sid, Worker* worker): sid(sid), worker(worker){}
		bool isnone() const{
			return sid.size() == 0 || worker == NULL;
		}
		//>>$ObjectA
	};
	struct Object2D: public Object{
		Object2D(std::string sid, Worker* worker): Object(sid, worker){}
		//>>$Object2D
	};
	struct Object3D: public Object{
		Object3D(std::string sid, Worker* worker): Object(sid, worker){}
		//>>$Object3D
	};
	struct Contour2D: public Object2D{
		Contour2D(std::string sid, Worker* worker): Object2D(sid, worker){}
		static Contour2D None(){
			return Contour2D("", NULL);
		}
		//>>$Contour2D
	};

	struct Grid2D: public Object2D{
		Grid2D(std::string sid, Worker* worker): Object2D(sid, worker){}
		static Grid2D None(){
			return Grid2D("", NULL);
		}
		//>>$Grid2D
	};

	struct Surface3D: public Object3D{
		Surface3D(std::string sid, Worker* worker): Object3D(sid, worker){}
		static Surface3D None(){
			return Surface3D("", NULL);
		}
		//>>$Surface3D
	};

	struct Grid3D: public Object3D{
		Grid3D(std::string sid, Worker* worker): Object3D(sid, worker){}
		static Grid3D None(){
			return Grid3D("", NULL);
		}
		//>>$Grid3D
	};

	//>>$Hybmesh
};

// ====================== Worker implementation
int Hybmesh::Worker::_to_int(const std::string& s){
	int ret;
	std::istringstream(s)>>ret;
	return ret;
}
double Hybmesh::Worker::_to_double(const std::string& s){
	double ret;
	std::istringstream(s)>>ret;
	return ret;
}
Hybmesh::Point2 Hybmesh::Worker::_to_point(const std::string& s){
	if (s == "None") return Hybmesh::Point2::None();
	double x, y; char _tmp;
	std::istringstream(s)>>_tmp>>x>>_tmp>>y;
	return Hybmesh::Point2(x, y);
}
std::vector<int> Hybmesh::Worker::_to_vecint(const std::string& s){
	std::vector<std::string> ss = __parse_vecstring(s);
	std::vector<int> ret(ss.size());
	for (size_t i=0; i<ret.size(); ++i){
		std::istringstream(ss[i])>>ret[i];
	}
	return ret;
}

//removes ending []; splits by ,; strips substrings from """, "'", " ";
std::vector<std::string> Hybmesh::Worker::__parse_vecstring(const std::string& s){
	std::vector<std::string> ret, ret2;
	if (s == "[]" || s == "None" || s.size() == 0) return ret;
	std::string::const_iterator it = s.begin();
	std::string::const_iterator itend = s.end() - 1;
	//remove []
	if (*it == '[') ++it;
	if (*itend == ']') --itend;
	//split
	std::string::const_iterator it0 = it;
	int bracket_level = 0;
	while (it <= itend){
		if (*it == '['){
			++bracket_level;
		} else if (*it == ']'){
			--bracket_level;
		} else if (bracket_level == 0){
			if (*it == ','){
				ret.push_back(std::string(it0, it));
				it0 = it + 1;
			} else if (it == itend){
				ret.push_back(std::string(it0, it+1));
			} else if (*it == ' '){
				it0 = it + 1;
			}
		}
		++it;
	}
	//strip
	for (size_t i=0; i<ret.size(); ++i){
		size_t i0 = ret[i].find_first_not_of("\'\" ");
		size_t i1 = ret[i].find_last_not_of("\'\" ");
		if (i0 <= i1) ret2.push_back(ret[i].substr(i0, i1 - i0 + 1));
	}
	return ret2;
}
template<class Obj> Obj Hybmesh::Worker::_to_object(const std::string& s){
	std::vector<std::string> ss = __parse_vecstring(s);
	if (ss.size() == 0 || ss[0] == "None"){
		return Obj::None();
	} else{
		return Obj(ss[0], this);
	}
}
template<class Obj> std::vector<Obj> Hybmesh::Worker::_to_vecobject(const std::string& s){
	std::vector<std::string> ss = __parse_vecstring(s);
	if (ss.size() == 0){
		return std::vector<Obj>();
	} else{
		std::vector<Obj> ret;
		for (size_t i=0; i<ss.size(); ++i){
			ret.push_back(_to_object<Obj>(ss[i]));
		}
		return ret;
	}
}
Hybmesh::Grid2D Hybmesh::Worker::_to_grid(const std::string& s){
	return _to_object<Hybmesh::Grid2D>(s);
}
Hybmesh::Contour2D Hybmesh::Worker::_to_cont(const std::string& s){
	return _to_object<Hybmesh::Contour2D>(s);
}
Hybmesh::Grid3D Hybmesh::Worker::_to_grid3(const std::string& s){
	return _to_object<Hybmesh::Grid3D>(s);
}
Hybmesh::Surface3D Hybmesh::Worker::_to_surface(const std::string& s){
	return _to_object<Hybmesh::Surface3D>(s);
}
std::vector<Hybmesh::Contour2D> Hybmesh::Worker::_to_veccont(const std::string& s){
	return _to_vecobject<Hybmesh::Contour2D>(s);
}
std::vector<Hybmesh::Grid2D> Hybmesh::Worker::_to_vecgrid(const std::string& s){
	return _to_vecobject<Hybmesh::Grid2D>(s);
}
std::vector<Hybmesh::Surface3D> Hybmesh::Worker::_to_vecsurface(const std::string& s){
	return _to_vecobject<Hybmesh::Surface3D>(s);
}
std::vector<Hybmesh::Grid3D> Hybmesh::Worker::_to_vecgrid3(const std::string& s){
	return _to_vecobject<Hybmesh::Grid3D>(s);
}

// cpp type -> string
std::string Hybmesh::Worker::_tos_bool(bool val){
	return val ? "True" : "False";
}
std::string Hybmesh::Worker::_tos_int(int val){
	std::ostringstream os;
	os<<val;
	return os.str();
}
std::string Hybmesh::Worker::_tos_double(double val){
	std::ostringstream os; os.precision(17);
	os<<val;
	return os.str();
}
std::string Hybmesh::Worker::_tos_string(const std::string& val){
	std::ostringstream os;
	os<<'"'<<val<<'"';
	return os.str();
}
std::string Hybmesh::Worker::_tos_point(Point2 val){
	if (val.isnone()) return "None";
	std::ostringstream os; os.precision(17);
	os<<'['<<val.x<<", "<<val.y<<']';
	return os.str();
}
std::string Hybmesh::Worker::_tos_point3(Point3 val){
	if (val.isnone()) return "None";
	std::ostringstream os; os.precision(17);
	os<<'['<<val.x<<", "<<val.y<<", "<<val.z<<']';
	return os.str();
}
std::string Hybmesh::Worker::_tos_vecbyte(const VecByte& comstr){
	std::string ret(comstr.begin(), comstr.end());
	while (ret.size() > 0 && ret.back() == '\0') ret.resize(ret.size() - 1);
	return ret;
}
std::string Hybmesh::Worker::_tos_vecint(const std::vector<int>& val){
	std::ostringstream os;
	os<<'[';
	for (size_t i=0; i<val.size(); ++i){
		os<<val[i]<<", ";
	}
	os<<']';
	return os.str();
}
std::string Hybmesh::Worker::_tos_vecdouble(const std::vector<double>& val){
	std::ostringstream os; os.precision(17);
	os<<'[';
	for (size_t i=0; i<val.size(); ++i){
		os<<val[i]<<", ";
	}
	os<<']';
	return os.str();
}
std::string Hybmesh::Worker::_tos_vecstring(const std::vector<std::string>& val){
	std::ostringstream os;
	os<<'[';
	for (size_t i=0; i<val.size(); ++i){
		os<<'"'<<val[i]<<"\", ";
	}
	os<<']';
	return os.str();
}
std::string Hybmesh::Worker::_tos_vecpoint(const std::vector<Point2>& val){
	std::ostringstream os; os.precision(17);
	os<<'[';
	for (size_t i=0; i<val.size(); ++i){
		os<<_tos_point(val[i])<<", ";
	}
	os<<']';
	return os.str();
}

template<class Obj>
std::string Hybmesh::Worker::_tos_object(const Obj& val){
	if (val.isnone()) return "None";
	std::ostringstream ss;
	ss<<'"'<<val.sid<<'"';
	return ss.str();
}
template<class Obj>
std::string Hybmesh::Worker::_tos_vecobject(const std::vector<Obj>& val){
	std::vector<std::string> rs(val.size());
	for (size_t i=0; i<val.size(); ++i) rs[i] = val[i].sid;
	return _tos_vecstring(rs);
}
std::vector<std::pair<int, double> >
Hybmesh::Worker::_to_vec_int_double_raw(const VecByte& val){
	const char* ps = &val[0];
	int sz = *(int*)ps;
	ps += sizeof(int);
	std::vector<std::pair<int, double> > ret(sz);
	for (int i=0; i<sz; ++i){
		ret[i].first = *(int*)ps;
		ps += sizeof(int);
		ret[i].second = *(double*)ps;
		ps += sizeof(double);
	}
	return ret;
};
std::vector<double> Hybmesh::Worker::_to_vecdouble_raw(const VecByte& val){
	const char* ps = &val[0];
	int sz = *(int*)ps;
	ps += sizeof(int);
	std::vector<double> ret(sz);
	for (int i=0; i<sz; ++i){
		ret[i] = *(double*)ps;
		ps += sizeof(double);
	}
	return ret;
};

std::vector<int> Hybmesh::Worker::_to_vecint_raw(const VecByte& val){
	const char* ps = &val[0];
	int sz = *(int*)ps;
	ps += sizeof(int);
	std::vector<int> ret(sz);
	for (int i=0; i<sz; ++i){
		ret[i] = *(int*)ps;
		ps += sizeof(int);
	}
	return ret;
}

//command interface
Hybmesh::VecByte Hybmesh::Worker::_apply_command(
		const std::string& func, const std::string& com){
	//func - name of the function
	//com - comma separated list of arguments. no ().
	_send_command(func, com);
	while (true){
		char sig = _wait_for_signal();
		switch (sig){
			//normal return
			case 'R': return _read_buffer();
			//callback
			case 'B': _apply_callback(_read_buffer()); break;
			//exception return
			case 'I': throw Hybmesh::EUserInterrupt();
			case 'E': throw Hybmesh::ERuntimeError(
					_tos_vecbyte(_read_buffer()).c_str());
			//something went wrong
			case '0': throw Hybmesh::ERuntimeError(
					"Server stopped working");
			default: throw Hybmesh::ERuntimeError(
					"Invalid client instruction");
		}
	};
}

void Hybmesh::Worker::_apply_callback(const VecByte& buf){
	const char* ps = &buf[0];
	double p1 = *(double*)(ps); ps+=sizeof(double);
	double p2 = *(double*)(ps); ps+=sizeof(double);
	int len1 = *(int*)(ps); ps+=sizeof(int);
	int len2 = *(int*)(ps); ps+=sizeof(int);
	const char* n1 = ps; ps+=len1;
	const char* n2 = ps; ps+=len2;
	std::string s1(n1, n1+len1);
	std::string s2(n2, n2+len2);
	int r = callback(s1, s2, p1, p2);
	char ret = (r == 0) ? 'G' : 'S';
	send_signal(ret);
}

void Hybmesh::Worker::_send_command(const std::string& func, const std::string& com){
	send_data(func.size(), &func[0]);
	send_data(com.size(), &com[0]);
	send_signal('C');
}

char Hybmesh::Worker::_wait_for_signal(){
	char ret;
	get_signal(&ret);
	return ret;
}

Hybmesh::VecByte Hybmesh::Worker::_read_buffer(){
	char* buf;
	int sz;
	get_data(&sz, &buf);
	Hybmesh::VecByte ret(buf, buf+sz);
	free(buf);
	return ret;
}
//=============== low level functionality
#include "unistd.h"
#ifdef WIN32

#include "windows.h"
#include "io.h"
#include "fcntl.h"

#define HMPLATFORM_READ _read
#define HMPLATFORM_WRITE _write
#define HMPLATFORM_CLOSE _close

#else

#include "sys/wait.h"
#define HMPLATFORM_READ read
#define HMPLATFORM_WRITE write
#define HMPLATFORM_CLOSE close

#endif

void Hybmesh::Worker::require_connection(const char* path){
#ifdef WIN32 // ====================== WINDOWS IMPLEMENTATION
	char cmd[1000];
	HANDLE sig_client2server[2];
	HANDLE sig_server2client[2];
	HANDLE data_client2server[2];
	HANDLE data_server2client[2];
	SECURITY_ATTRIBUTES sa;
	PROCESS_INFORMATION pi;
	STARTUPINFO si;

	/* create pipes */
	sa.nLength = sizeof(sa);
	sa.bInheritHandle = TRUE;
	sa.lpSecurityDescriptor = NULL;
	if (!CreatePipe(&sig_client2server[0], &sig_client2server[1], &sa, 0) || 
	    !CreatePipe(&sig_server2client[0], &sig_server2client[1], &sa, 0) ||
	    !CreatePipe(&data_client2server[0], &data_client2server[1], &sa, 0) ||
	    !CreatePipe(&data_server2client[0], &data_server2client[1], &sa, 0)){
		return 0;
	}
	/* do not inherit client handles */
	SetHandleInformation(sig_client2server[1], HANDLE_FLAG_INHERIT, 0);
	SetHandleInformation(sig_server2client[0], HANDLE_FLAG_INHERIT, 0);
	SetHandleInformation(data_client2server[1], HANDLE_FLAG_INHERIT, 0);
	SetHandleInformation(data_server2client[0], HANDLE_FLAG_INHERIT, 0);
	SetHandleInformation(sig_client2server[0], HANDLE_FLAG_INHERIT, 1);
	SetHandleInformation(sig_server2client[1], HANDLE_FLAG_INHERIT, 1);
	SetHandleInformation(data_client2server[0], HANDLE_FLAG_INHERIT, 1);
	SetHandleInformation(data_server2client[1], HANDLE_FLAG_INHERIT, 1);

	/* forking */
	ZeroMemory(&si, sizeof(si));
	ZeroMemory(&pi, sizeof(pi));
	si.cb = sizeof(si);

	intptr_t h1 = (intptr_t)sig_client2server[0];
	intptr_t h2 = (intptr_t)sig_server2client[1];
	intptr_t h3 = (intptr_t)data_client2server[0];
	intptr_t h4 = (intptr_t)data_server2client[1];

	sprintf(cmd, "%s %lld %lld %lld %lld", path, h1, h2, h3, h4);

	if(!CreateProcess( NULL,/* No module name (use command line) */
		cmd,            /* Command line*/
		NULL,           /* Process handle not inheritable */
		NULL,           /* Thread handle not inheritable */
		TRUE,           /* Set handle inheritance to TRUE */
		0,              /* No creation flags */
		NULL,           /* Use parent's environment block */
		NULL,           /* Use parent's starting directory */
		&si,            /* Pointer to STARTUPINFO structure */
		&pi )           /* Pointer to PROCESS_INFORMATION structure */
	){
		std::ostringstream errs;
		errs<<"CreateProcess failed with error "<<GetLastError();
		throw Hybmesh::ERuntimeError(errs.str().c_str());
	};

	/* fill con */
	this->childid = (intptr_t)pi.hProcess;
	this->sig_read = _open_osfhandle((intptr_t)sig_server2client[0], _O_APPEND|_O_RDONLY);
	this->sig_write = _open_osfhandle((intptr_t)sig_client2server[1], _O_APPEND);
	this->data_read = _open_osfhandle((intptr_t)data_server2client[0], _O_APPEND|_O_RDONLY);
	this->data_write = _open_osfhandle((intptr_t)data_client2server[1], _O_APPEND);

	/* close unused handles */
	CloseHandle(pi.hThread);
	CloseHandle(sig_server2client[1]);
	CloseHandle(sig_client2server[0]);
	CloseHandle(data_server2client[1]);
	CloseHandle(data_client2server[0]);

#else // ==================== POSIX IMPLEMENTATION
	char s0[16], s1[16], s2[16], s3[16];
	int sig_client2server[2],  sig_server2client[2],
	    data_client2server[2], data_server2client[2];
	int childid;

	/* create pipes */
	if (pipe(sig_client2server)<0 || pipe(sig_server2client)<0 ||
	    pipe(data_client2server)<0 || pipe(data_server2client)<0){
		throw Hybmesh::ERuntimeError("failed to establish a pipe");
	}

	/* forking */
	childid = fork();
	if (childid < 0) throw Hybmesh::ERuntimeError("failed to fork");
	else if (childid == 0){
		/* child process: server */
		sprintf(s0, "%d", sig_client2server[0]);
		sprintf(s1, "%d", sig_server2client[1]);
		sprintf(s2, "%d", data_client2server[0]);
		sprintf(s3, "%d", data_server2client[1]);
		close(sig_client2server[1]);
		close(sig_server2client[0]);
		close(data_client2server[1]);
		close(data_server2client[0]);
		int err = execl(path, path, "-px", s0, s1, s2, s3, NULL);
		if (err == -1) throw ERuntimeError(
			"failed to launch hybmesh server application");
		return;
	} else {
		/* parent process: client */
		close(sig_client2server[0]);
		close(sig_server2client[1]);
		close(data_client2server[0]);
		close(data_server2client[1]);
		this->childid = childid;
		this->sig_read = sig_server2client[0];
		this->sig_write = sig_client2server[1];
		this->data_read = data_server2client[0];
		this->data_write = data_client2server[1];
	}

#endif //POSIX
}


void Hybmesh::Worker::get_signal(char* sig){
	*sig = '0';
	HMPLATFORM_READ(sig_read, sig, 1);
}
void Hybmesh::Worker::get_data(int* sz, char** data){
	char intbuf[4];
	HMPLATFORM_READ(data_read, intbuf, 4);
	*sz = *(int*)intbuf;
	*data = (char*)malloc(*sz);
	HMPLATFORM_READ(data_read, *data, *sz);
}
void Hybmesh::Worker::send_signal(char sig){
	HMPLATFORM_WRITE(sig_write, &sig, 1);
}
void Hybmesh::Worker::send_data(int sz, const char* data){
	HMPLATFORM_WRITE(data_write, (char*)(&sz), 4);
	HMPLATFORM_WRITE(data_write, data, sz);
}
void Hybmesh::Worker::break_connection(){
	send_signal('Q');
#ifdef WIN32
	WaitForSingleObject((HANDLE)childid, INFINITE);
#else
	waitpid(childid, NULL, 0);
#endif
	HMPLATFORM_CLOSE(sig_read);
	HMPLATFORM_CLOSE(sig_write);
	HMPLATFORM_CLOSE(data_read);
	HMPLATFORM_CLOSE(data_write);
}


#endif
