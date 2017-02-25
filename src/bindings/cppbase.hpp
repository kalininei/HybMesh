#ifndef HYBMESH_CPP_FRONTEND_HPP
#define HYBMESH_CPP_FRONTEND_HPP

#include <vector>
#include <string>
#include <limits>

class Hybmesh{
public:
	struct Point2;
	struct Point3;
	struct Contour2D;
	struct Grid2D;
	struct Surface3D;
	struct Grid3D;
	typedef std::vector<char> VecByte;
private:

	struct Worker{
		VecByte _apply_command(const std::string& func, const std::string& com);
		// string -> cpp type
		int _to_int(const std::string&);
		double _to_double(const std::string&);
		Point2 _to_point(const std::string&);
		Grid2D _to_grid(const std::string&);
		Contour2D _to_cont(const std::string&);
		Grid3D _to_grid3(const std::string&);
		std::vector<int> _to_vecint(const std::string&);
		std::vector<Contour2D> _to_veccont(const std::string&);
		std::vector<Grid2D> _to_vecgrid(const std::string&);
		std::vector<Surface3D> _to_vecsurface(const std::string&);
		std::vector<Grid3D> _to_vecgrid3(const std::string&);
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
		template<class Obj>
		std::string _tos_vecobject(const std::vector<Obj>&);
		// vecbyte -> cpp type
		std::vector<std::pair<int, double> > _to_map_int_double_raw(const VecByte&);
		std::vector<double> _to_vecdouble_raw(const VecByte&);
	};

	Worker _worker;
	Worker* worker;

	// copy is forbidden
	Hybmesh(const Hybmesh&);
public:
	struct Point2{
		double x, y;
		Point2(double x=0, double y=0): x(x), y(y){}

		static Point2 None(){
			return Point2(std::numeric_limits<double>::max(),
					std::numeric_limits<double>::max());
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
	};

	struct Object{
		std::string sid;
		Worker* worker;
		Object(std::string sid, Worker* worker): sid(sid), worker(worker){}
	};
	struct Object2D: public Object{
		Object2D(std::string sid, Worker* worker): Object(sid, worker){}
	};
	struct Object3D: public Object{
		Object3D(std::string sid, Worker* worker): Object(sid, worker){}
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

#endif
