using System;
using System.Collections.Generic;

public class Hybmesh{
//private:
	private class Worker{
		public byte[] _apply_command( String func,  String com){
			throw new NotImplementedException();
		}
		// string -> cpp type
		public int _to_int(String str){
			throw new NotImplementedException();
		}
		public double _to_double(String str){
			throw new NotImplementedException();
		}
		public Point2 _to_point(String str){
			throw new NotImplementedException();
		}
		public Grid2D _to_grid(String str){
			throw new NotImplementedException();
		}
		public Contour2D _to_cont(String str){
			throw new NotImplementedException();
		}
		public Grid3D _to_grid3(String str){
			throw new NotImplementedException();
		}
		public int[] _to_vecint(String str){
			throw new NotImplementedException();
		}
		public Contour2D[] _to_veccont(String str){
			throw new NotImplementedException();
		}
		public Grid2D[] _to_vecgrid(String str){
			throw new NotImplementedException();
		}
		public Surface3D[] _to_vecsurface(String str){
			throw new NotImplementedException();
		}
		public Grid3D[] _to_vecgrid3(String str){
			throw new NotImplementedException();
		}
		// cpp type -> string
		public String _tos_bool(bool val){
			throw new NotImplementedException();
		}
		public String _tos_int(int val){
			throw new NotImplementedException();
		}
		public String _tos_double(doub val){
			throw new NotImplementedException();
		}
		public String _tos_point(Point2 val){
			throw new NotImplementedException();
		}
		public String _tos_point3(Point3 val){
			throw new NotImplementedException();
		}
		public String _tos_vecbyte(byte[] val){
			throw new NotImplementedException();
		}
		public String _tos_vecint(int[] val){
			throw new NotImplementedException();
		}
		public String _tos_vecdouble(double[] val){
			throw new NotImplementedException();
		}
		public String _tos_vecstring(String[] val){
			throw new NotImplementedException();
		}
		public String _tos_vecpoint(Point2[] val){
			throw new NotImplementedException();
		}
		public String _tos_vecobject<Obj>(Obj[] val){
			throw new NotImplementedException();
		}
		// vecbyte -> cpp type
		public KeyValuePair<int, double> _to_map_int_double_raw(byte[] val){
			throw new NotImplementedException();
		}
		public double[] _to_vecdouble_raw(byte[] val){
			throw new NotImplementedException();
		}
	};

	Worker worker;

//public:
	public class Point2{
		public double x, y;
		public Point2(double x=0, double y=0){
			this.x=x;
			this.y=y;
		}
	};

	public class Point3{
		public double x, y, z;
		public Point3(double x=0, double y=0, double z=0){
			this.x=x;
			this.y=y;
			this.z=z;
		}
	};

	public class Object{
		public String sid;
		public Worker worker;
		public Object(String sid, Worker worker){
			this.sid=sid;
			this.worker=worker;
		}
	};
	public class Object2D: Object{
		public Object2D(String sid, Worker worker): base(sid, worker){}
	};
	public class Object3D: Object{
		public Object3D(String sid, Worker worker): base(sid, worker){}
	};

	public class Contour2D: Object2D{
		public Contour2D(String sid, Worker worker): base(sid, worker){}
		//>>$Contour2D
	};

	public class Grid2D: Object2D{
		public Grid2D(String sid, Worker worker): base(sid, worker){}
		//>>$Grid2D
	};

	public class Surface3D: Object3D{
		public Surface3D(String sid, Worker worker): base(sid, worker){}
		//>>$Surface3D
	};

	public class Grid3D: Object3D{
		public Grid3D(String sid, Worker worker): base(sid, worker){}
		//>>$Grid3D
	};

	//>>$Hybmesh
};
