import java.util.HashMap;
import java.util.AbstractMap;

class Hybmesh{
//private:
	private class Worker{
		public byte[] _apply_command(String func, String com){
			throw new RuntimeException("Not Implemented");
		}
		// string -> cpp type
		public int _to_int(String str){
			throw new RuntimeException("Not implemented");
		}
		public double _to_double(String str){
			throw new RuntimeException("Not implemented");
		}
		public Point2 _to_point(String str){
			throw new RuntimeException("Not implemented");
		}
		public Grid2D _to_grid(String str){
			throw new RuntimeException("Not implemented");
		}
		public Contour2D _to_cont(String str){
			throw new RuntimeException("Not implemented");
		}
		public Grid3D _to_grid3(String str){
			throw new RuntimeException("Not implemented");
		}
		public int[] _to_vecint(String str){
			throw new RuntimeException("Not implemented");
		}
		public Contour2D[] _to_veccont(String str){
			throw new RuntimeException("Not implemented");
		}
		public Grid2D[] _to_vecgrid(String str){
			throw new RuntimeException("Not implemented");
		}
		public Surface3D[] _to_vecsurface(String str){
			throw new RuntimeException("Not implemented");
		}
		public Grid3D[] _to_vecgrid3(String str){
			throw new RuntimeException("Not implemented");
		}
		// cpp type -> string
		public String _tos_bool(boolean val){
			throw new RuntimeException("Not implemented");
		}
		public String _tos_int(int val){
			throw new RuntimeException("Not implemented");
		}
		public String _tos_double(double val){
			throw new RuntimeException("Not implemented");
		}
		public String _tos_point(Point2 val){
			throw new RuntimeException("Not implemented");
		}
		public String _tos_point3(Point3 val){
			throw new RuntimeException("Not implemented");
		}
		public String _tos_vecbyte(byte[] val){
			throw new RuntimeException("Not implemented");
		}
		public String _tos_vecint(int[] val){
			throw new RuntimeException("Not implemented");
		}
		public String _tos_vecdouble(double[] val){
			throw new RuntimeException("Not implemented");
		}
		public String _tos_vecstring(String[] val){
			throw new RuntimeException("Not implemented");
		}
		public String _tos_vecpoint(Point2[] val){
			throw new RuntimeException("Not implemented");
		}
		public <Obj> String _tos_vecobject(Obj[] val){
			throw new RuntimeException("Not implemented");
		}
		// vecbyte -> cpp type
		public AbstractMap.SimpleEntry<Integer, Double>[] _to_map_int_double_raw(byte[] val){
			throw new RuntimeException("Not implemented");
		}
		public double[] _to_vecdouble_raw(byte[] val){
			throw new RuntimeException("Not implemented");
		}
	};

	Worker worker;
//public:
	public class Point2{
		public double x, y;
		public Point2(double x, double y){
			this.x = x;
			this.y = y;
		};
	};
	public class Point3{
		public double x, y, z;
		public Point3(double x, double y, double z){
			this.x = x;
			this.y = y;
			this.z = z;
		};
	};
	public class Object{
		public String sid;
		public Worker worker;
		public Object(String sid_, Worker worker_){
			sid = sid_;
			worker = worker_;
		}
	};
	public class Object2D extends Object{
		public Object2D(String sid, Worker worker){
			super(sid, worker);
		}
	};
	public class Object3D extends Object{
		public Object3D(String sid, Worker worker){
			super(sid, worker);
		}
	};

	public class Contour2D extends Object2D{
		public Contour2D(String sid, Worker worker){
			super(sid, worker);
		};
		//>>$Contour2D
	};

	public class Grid2D extends Object2D{
		public Grid2D(String sid, Worker worker){
			super(sid, worker);
		};
		//>>$Grid2D
	};

	public class Surface3D extends Object3D{
		public Surface3D(String sid, Worker worker){
			super(sid, worker);
		}
		//>>$Surface3D
	};

	public class Grid3D extends Object3D{
		public Grid3D(String sid, Worker worker){
			super(sid, worker);
		}
		//>>$Grid3D
	};

	//>>$Hybmesh
};
