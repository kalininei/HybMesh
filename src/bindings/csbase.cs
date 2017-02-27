using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;

internal class Worker: IDisposable{
	// ==== communication library
	private int connection = -1;
	[DllImport(hmpipespath)]
	private static extern int require_connection(
		[MarshalAs(UnmanagedType.LPStr)] string server_path);
	[DllImport(hmpipespath)]
	private static extern byte get_signal(int con);
	[DllImport(hmpipespath)]
	private static extern int get_data(int con, out IntPtr data);
	[DllImport(hmpipespath)]
	private static extern void send_signal(int con, byte sig);
	[DllImport(hmpipespath)]
	private static extern void send_data(int con, int sz, byte[] data);
	[DllImport(hmpipespath)]
	private static extern void break_connection(int con);
	[DllImport(hmpipespath)]
	private static extern void free_char_array(IntPtr s);
	//server communication
	private void _send_command(string func, string com){
		byte[] dtfunc = Encoding.UTF8.GetBytes(func);
		byte[] dtcom = Encoding.UTF8.GetBytes(com);
		send_data(connection, dtfunc.Length, dtfunc);
		send_data(connection, dtcom.Length, dtcom);
		send_signal(connection, (byte)'C');
	}
	private byte _wait_for_signal(){
		return get_signal(connection);
	}
	private byte[] _read_buffer(){
		IntPtr buf;
		int sz = get_data(connection, out buf);
		byte[] ret = new byte[sz];
		Marshal.Copy(buf, ret, 0, sz);
		free_char_array(buf);
		return ret;
	}
	//command interface
	public byte[] _apply_command(string func,  string com){
		_send_command(func, com);
		while (true){
			byte sig = _wait_for_signal();
			switch (sig){
				//normal return
				case (byte)'R':
					return _read_buffer();
				//callback
				case (byte)'B': 
					_apply_callback(_read_buffer());
					break;
				//exception return
				case (byte)'I': 
					throw new Hybmesh.EUserInterrupt();
				case (byte)'E':
					throw new Hybmesh.ERuntimeError(
						_tos_vecbyte(_read_buffer()));
				//something went wrong
				case (byte)'0': 
					throw new Hybmesh.ERuntimeError(
						"Server stopped working");
				default:
					throw new Hybmesh.ERuntimeError(
						"Invalid client instruction");
			}
		};
	}
	public void  _apply_callback(byte[] buf){
		using (var mstream = new MemoryStream(buf))
		using (var reader = new BinaryReader(mstream)){
			double p1 = reader.ReadDouble();
			double p2 = reader.ReadDouble();
			int len1 = reader.ReadInt32();
			int len2 = reader.ReadInt32();
			char[] n1 = reader.ReadChars(len1);
			char[] n2 = reader.ReadChars(len2);
			string s1 = new string(n1);
			string s2 = new string(n2);
			int r = callback(s1, s2, p1, p2);
			byte ret = (r == 0) ? (byte)'G' : (byte)'S';
			send_signal(connection, ret);
		}
	}
	//string -> cpp type
	private string[] __parse_vecstring(string str){
		string s2 = str.Substring(1, str.Length-2);
		return str.Split(",");
	}
	private Obj[] __parse_vecobject<Obj>(string str){
		Obj[] ret = new Obj[]{};
		string[] s2 = __parse_vecstring(str);
		foreach(var it in s2) ret.Add(new Obj(it, this));
		return ret;
	}
	public int _to_int(string str){
		return Int32.Parse(str);
	}
	public double _to_double(string str){
		return double.Parse(str);
	}
	public Point2 _to_point(string str){
		if (str == "None") return null;
		string[] s3 = __parse_vecstring(str);
		return new Point(double.Parse(s3[0]), double.Parse(s3[1]));
	}
	public Grid2D _to_grid(string str){
		if (str == "None") return null;
		else if (str[0] == '[') return _to_vecgrid(str)[0];
		else return new Grid2D(str, this);
	}
	public Contour2D _to_cont(string str){
		if (str == "None") return null;
		else if (str[0] == '[') return _to_veccont(str)[0];
		else return new Contour2D(str, this);
	}
	public Grid3D _to_grid3(string str){
		if (str == "None") return null;
		else if (str[0] == '[') return _to_vecgrid3(str)[0];
		else return new Grid3D(str, this);
	}
	public int[] _to_vecint(string str){
		string[] ss = __parse_vecstring(str);
		int[] ret = new int[ss.Length];
		for (int i=0; i<ret.Length; ++i){
			ret[i] = _to_int(ss[i]);
		}
		return ret;
	}
	public Contour2D[] _to_veccont(string str){
		if (str[0] != '[') return new Contour2D[]{ _to_cont(str) };
		else return __parse_vecobj<Contour2D>(str);
	}
	public Grid2D[] _to_vecgrid(string str){
		if (str[0] != '[') return new Grid2D[]{ _to_grid(str) };
		else return __parse_vecobj<Grid2D>(str);
	}
	public Surface3D[] _to_vecsurface(string str){
		if (str[0] != '[') return new Surface3D[]{ _to_surface(str) };
		else return __parse_vecobj<Surface3D>(str);
	}
	public Grid3D[] _to_vecgrid3(string str){
		if (str[0] != '[') return new Grid3D[]{ _to_grid3(str) };
		else return __parse_vecobj<Grid3D>(str);
	}
	//cpp type -> string
	public string _tos_bool(bool val){
		return val ? "True" : "False";
	}
	public string _tos_int(int val){
		return val.Tostring();
	}
	public string _tos_double(double val){
		return val.Tostring("G17", CultureInfo.InvariantCulture);
	}
	public string _tos_point(Point2 val){
		if (val == null) return "None";
		else return '[' + val.x.Tostring("G17", CultureInfo.InvariantCulture) +
			',' + val.y.Tostring("G17", CultureInfo.InvariantCulture) +
			']';
	}
	public string _tos_point3(Point3 val){
		if (val == null) return "None";
		else return '[' + val.x.Tostring("G17", CultureInfo.InvariantCulture) +
			',' + val.y.Tostring("G17", CultureInfo.InvariantCulture) +
			',' + val.z.Tostring("G17", CultureInfo.InvariantCulture) +
			']';
	}
	public string _tos_vecbyte(byte[] val){
		return Encoding.UTF8.Getstring(val);
	}
	public string _tos_vecint(int[] val){
		return '[' + string.Join(", ", val.Select(
				p=>p.Tostring().ToArray())) + ']';
	}
	public string _tos_vecdouble(double[] val){
		return '[' + string.Join(", ", val.Select(
			p=>p.Tostring("G17", CultureInfo.InvariantCulture)).ToArray())
			+ ']';
	}
	public string _tos_vecstring(string[] val){
		return '[' + string.Join(", ", val.Select(
			p=>('"' + p + '"')).ToArray()) + ']';
	}
	public string _tos_vecpoint(Point2[] val){
		return '[' + string.Join(", ", val.Select(
			p=>this._tos_point(p)).ToArray()) + ']';
	}
	public string _tos_object<Obj>(Obj val){
		return '"' + val.sid + '"';
	}
	public string _tos_vecobject<Obj>(Obj[] val){
		return '[' + string.Join(", ", val.Select(
			p=>this._tos_object<Obj>(p)).ToArray()) + ']';
	}
	//vecbyte -> cpp type
	public KeyValuePair<int, double> _to_vec_int_double_raw(byte[] val){
		using (var mstream = new MemoryStream(val))
		using (var reader = new BinaryReader(mstream)){
			int len = reader.ReadInt32();
			var ret = new KeyValuePair<int, double>[]{};
			for (int i=0; i<len; ++i){
				int k = reader.ReadInt32();
				double v = reader.ReadDouble();
				ret.Add(new KeyValuePair<int, double>(k, v));
			}
			return ret;
		}
	}
	public double[] _to_vecdouble_raw(byte[] val){
		using (var mstream = new MemoryStream(val))
		using (var reader = new BinaryReader(mstream)){
			int len = reader.ReadInt32();
			var ret = new double[len];
			for (int i=0; i<len; ++i){
				ret[i] = reader.ReadDouble();
			}
			return ret;
		}
	}
	//callback
	public Hybmesh.DCallback callback;
	//constructor/destructor
	public Worker(){
		connection = require_connection(hmserverpath);
		callback = (string s1, string s2, double p1, double p2) => 1;
	}
	protected virtual void Dispose(bool disposing){
		if (connection != -1){
			break_connection(connection);
			connection = -1;
		}
		if (disposing){ /* release other disposable objects */ }
	}
	~Worker(){
		Dispose(false);
	}
	public void Dispose(){
		Dispose(true);
		GC.SuppressFinalize(this);
	}
};

public class Hybmesh: IDisposable{
	private Worker worker;
//public:
	public static string hybmesh_exec_path =
		"/home/ek/Hybmesh/build/bin/hybmesh";  //>>$EXEPATH
	public static string hybmesh_lib_path =
		"/home/ek/Hybmesh/build/lib/core_hybmesh_connection_cs.so";  //>>$LIBPATH

	[Serializable]
	public class EUserInterrupt: Exception{
		public EUserInterrupt(): base ("Interrupted by user");
	}
	[Serializable]
	public class ERuntimeError: Exception{
		public ERuntimeError(string message):
			base("Hybmesh runtime error:\n" + message){}
	};
	public delegate int DCallback(string s1, string s2, double p1, double p2);
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

	//=================== construction/destruction
	public Hybmesh(){
		worker = new Worker(hmserverpath);
	}
	public void Dispose(){
		worker.Dispose();
		GC.SuppressFinalize(this);
	}
	//==================== callback
	void AssignCallback(DCallback cb){
		worker.callback = cb;
	}

	// =================== Geometrical objects
	public class Object{
		public string sid;
		public Worker worker;
		public Object(string sid, Worker worker){
			this.sid=sid;
			this.worker=worker;
		}
	};
	public class Object2D: Object{
		public Object2D(string sid, Worker worker): base(sid, worker){}
	};
	public class Object3D: Object{
		public Object3D(string sid, Worker worker): base(sid, worker){}
	};

	public class Contour2D: Object2D{
		public Contour2D(string sid, Worker worker): base(sid, worker){}
		//>>$Contour2D
	};

	public class Grid2D: Object2D{
		public Grid2D(string sid, Worker worker): base(sid, worker){}
		//>>$Grid2D
	};

	public class Surface3D: Object3D{
		public Surface3D(string sid, Worker worker): base(sid, worker){}
		//>>$Surface3D
	};

	public class Grid3D: Object3D{
		public Grid3D(string sid, Worker worker): base(sid, worker){}
		//>>$Grid3D
	};

	//>>$Hybmesh
};
