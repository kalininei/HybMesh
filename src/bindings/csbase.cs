using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using System.Globalization;
using System.Text;
using System.IO;
using System.Linq;

/// <summary>
/// Provides interface for pipe communication with hybmesh executable
/// </summary>
public class Hybmesh: IDisposable{
	//Worker class and instance should be private.
	//But if so Object.worker instance can not be protected (visible for Grid2D, etc).
	
	public class Worker{
		// Use LoadLibrary/dlopen for core_hmconnection_cs  to be able to set its path on runtime
		[DllImport("libdl.so")]
		protected static extern IntPtr dlopen(string filename, int flags);
		[DllImport("libdl.so")]
		protected static extern IntPtr dlsym(IntPtr handle, string symbol);
		[DllImport("kernel32.dll")]
		protected static extern IntPtr LoadLibrary(string filename);
		[DllImport("kernel32.dll")]
		protected static extern IntPtr GetProcAddress(IntPtr hModule, string procname);

		public delegate int DRequireConnection(string s);
		public delegate byte DGetSignal(int con);
		public delegate int DGetData(int con, out IntPtr data);
		public delegate void DSendSignal(int con, byte sig);
		public delegate void DSendData(int con, int sz, byte[] data);
		public delegate void DBreakConnection(int con);
		public delegate void DFreeCharArray(IntPtr s);

		public static DRequireConnection require_connection;
		public static DGetSignal get_signal;
		public static DGetData get_data;
		public static DSendSignal send_signal;
		public static DSendData send_data;
		public static DBreakConnection break_connection;
		public static DFreeCharArray free_char_array;


		static Worker(){
			int p = (int) Environment.OSVersion.Platform;
			bool IsLinux = (p == 4) || (p == 6) || (p == 128);

			// 1. library name
			string libp;
			if (IsLinux){
				libp = Hybmesh.hybmesh_lib_path + "/libcore_hmconnection_cs.so";
			} else {
				libp = Hybmesh.hybmesh_lib_path + "/core_hmconnection_cs.dll";
			}

			//2. load library
			IntPtr moduleHandle;
			if (IsLinux){
				moduleHandle = dlopen(libp, 2);
			} else {
				moduleHandle = LoadLibrary(libp);
			}
			if (moduleHandle == IntPtr.Zero){
				throw new Exception("Failed to load library at " + libp);
			}

			// 3. load functions pointers
			IntPtr freq, fgetsig, fsendsig, fgetdata, fsenddata, fbreak, ffree;
			if (IsLinux){
				freq = dlsym(moduleHandle, "require_connection");
				fgetsig = dlsym(moduleHandle, "get_signal");
				fgetdata = dlsym(moduleHandle, "get_data");
				fsendsig = dlsym(moduleHandle, "send_signal");
				fsenddata = dlsym(moduleHandle, "send_data");
				fbreak = dlsym(moduleHandle, "break_connection");
				ffree = dlsym(moduleHandle, "free_char_array");
			} else {
				freq = GetProcAddress(moduleHandle, "require_connection");
				fgetsig = GetProcAddress(moduleHandle, "get_signal");
				fgetdata = GetProcAddress(moduleHandle, "get_data");
				fsendsig = GetProcAddress(moduleHandle, "send_signal");
				fsenddata = GetProcAddress(moduleHandle, "send_data");
				fbreak = GetProcAddress(moduleHandle, "break_connection");
				ffree = GetProcAddress(moduleHandle, "free_char_array");
			}

			//4. Assign delegates
			require_connection = Marshal.GetDelegateForFunctionPointer(
				freq, typeof(DRequireConnection)) as DRequireConnection;
			get_signal = Marshal.GetDelegateForFunctionPointer(
				fgetsig, typeof(DGetSignal)) as DGetSignal;
			get_data = Marshal.GetDelegateForFunctionPointer(
				fgetdata, typeof(DGetData)) as DGetData;
			send_signal = Marshal.GetDelegateForFunctionPointer(
				fsendsig, typeof(DSendSignal)) as DSendSignal;
			send_data = Marshal.GetDelegateForFunctionPointer(
				fsenddata, typeof(DSendData)) as DSendData;
			break_connection = Marshal.GetDelegateForFunctionPointer(
				fbreak, typeof(DBreakConnection)) as DBreakConnection;
			free_char_array = Marshal.GetDelegateForFunctionPointer(
				ffree, typeof(DFreeCharArray)) as DFreeCharArray;
		}
	
		// ==== communication library
		private int connection = -1;
		public void free(){
			if (connection != -1){
				break_connection(connection);
				connection = -1;
			}
		}
		//server communication
		private void _send_command(string func, string com){
			byte[] dtfunc = Encoding.UTF8.GetBytes(func);
			byte[] dtcom = Encoding.UTF8.GetBytes(com);
			send_signal(connection, (byte)'C');
			send_data(connection, dtfunc.Length, dtfunc);
			send_data(connection, dtcom.Length, dtcom);
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
		private void  _apply_callback(byte[] buf){
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
		//removes ending []; splits by ,; strips substrings from """, "'", " ";
		private string[] __parse_vecstring(string str){
			string s = str.Trim();
			if (s[0]=='[' && s[s.Length-1]==']') s = s.Substring(1, s.Length-2);
			if (s == "[]" || s == "None" || s.Length == 0) return new string[]{};
			//split by ',' keeping sublists [x, y] entries untouched
			List<string> ret1 = new List<string>();
			char[] ssymb = new char[]{'\'', '"', ' '};
			int bracket_level=0;
			foreach(var it in s.Split(',').Select(x=>x.Trim(ssymb))){
				if (it.Length==0) continue;
				else if (bracket_level == 0) ret1.Add(it);
				else ret1[ret1.Count()-1] = ret1[ret1.Count()-1]+','+it;
				if (it[0] == '[') ++bracket_level;
				if (it[it.Length-1] == ']') --bracket_level;
			}
			return ret1.ToArray();
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
			return new Point2(double.Parse(s3[0]), double.Parse(s3[1]));
		}
		private Obj _to_object<Obj>(string str) where Obj: Object{
			string[] ss = __parse_vecstring(str);
			if (ss.Length == 0 || ss[0] == "None"){
				return null;
			} else{
				return (Obj)Activator.CreateInstance(
						typeof(Obj), ss[0], this);
			}
		}
		private Obj[] _to_vecobj<Obj>(string str) where Obj: Object{
			string[] ss = __parse_vecstring(str);
			Obj[] ret = new Obj[ss.Length];
			for (int i=0; i<ss.Length; ++i){
				ret[i] = _to_object<Obj>(ss[i]);
			}
			return ret;
		}
		public Grid2D _to_grid(string str){
			return _to_object<Grid2D>(str);
		}
		public Contour2D _to_cont(string str){
			return _to_object<Contour2D>(str);
		}
		public Surface3D _to_surface(string str){
			return _to_object<Surface3D>(str);
		}
		public Grid3D _to_grid3(string str){
			return _to_object<Grid3D>(str);
		}
		public Contour2D[] _to_veccont(string str){
			return _to_vecobj<Contour2D>(str);
		}
		public Grid2D[] _to_vecgrid(string str){
			return _to_vecobj<Grid2D>(str);
		}
		public Surface3D[] _to_vecsurface(string str){
			return _to_vecobj<Surface3D>(str);
		}
		public Grid3D[] _to_vecgrid3(string str){
			return _to_vecobj<Grid3D>(str);
		}
		public int[] _to_vecint(string str){
			string[] ss = __parse_vecstring(str);
			int[] ret = new int[ss.Length];
			for (int i=0; i<ret.Length; ++i){
				ret[i] = _to_int(ss[i]);
			}
			return ret;
		}
		//cpp type -> string
		public string _tos_bool(bool val){
			return val ? "True" : "False";
		}
		public string _tos_int(int val){
			return val.ToString();
		}
		public string _tos_double(double val){
			return val.ToString("G16", CultureInfo.InvariantCulture);
		}
		public string _tos_point(Hybmesh.Point2 val){
			if (val == null) return "None";
			else return '[' + val.x.ToString("G16", CultureInfo.InvariantCulture) +
				',' + val.y.ToString("G16", CultureInfo.InvariantCulture) +
				']';
		}
		public string _tos_point3(Hybmesh.Point3 val){
			if (val == null) return "None";
			else return '[' + val.x.ToString("G16", CultureInfo.InvariantCulture) +
				',' + val.y.ToString("G16", CultureInfo.InvariantCulture) +
				',' + val.z.ToString("G16", CultureInfo.InvariantCulture) +
				']';
		}
		public string _tos_string(string val){
			if (val == null) return "None";
			return '"' + val + '"';
		}
		public string _tos_vecbyte(byte[] val){
			if (val == null) return "None";
			return Encoding.UTF8.GetString(val);
		}
		public string _tos_vecint(int[] val){
			if (val == null) return "None";
			return '[' + string.Join(", ", val.Select(
					p=>p.ToString()).ToArray()) + ']';
		}
		public string _tos_vecdouble(double[] val){
			if (val == null) return "None";
			return '[' + string.Join(", ", val.Select(
				p=>p.ToString("G16", CultureInfo.InvariantCulture)).ToArray())
				+ ']';
		}
		public string _tos_vecstring(string[] val){
			if (val == null) return "None";
			return '[' + string.Join(", ", val.Select(
				p=>('"' + p + '"')).ToArray()) + ']';
		}
		public string _tos_vecpoint(Hybmesh.Point2[] val){
			if (val == null) return "None";
			return '[' + string.Join(", ", val.Select(
				p=>this._tos_point(p)).ToArray()) + ']';
		}
		public string _tos_object(Object val){
			if (val == null) return "None";
			else return '"' + val.sid + '"';
		}
		public string _tos_vecobject(Object[] val){
			if (val == null) return "None";
			return '[' + string.Join(", ", val.Select(
				p=>this._tos_object(p)).ToArray()) + ']';
		}
		//vecbyte -> cpp type
		public KeyValuePair<int, double>[] _to_vec_int_double_raw(byte[] val){
			using (var mstream = new MemoryStream(val))
			using (var reader = new BinaryReader(mstream)){
				int len = reader.ReadInt32();
				var ret = new KeyValuePair<int, double>[len];
				for (int i=0; i<len; ++i){
					int k = reader.ReadInt32();
					double v = reader.ReadDouble();
					ret[i] = new KeyValuePair<int, double>(k, v);
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
		public int[] _to_vecint_raw(byte[] val){
			using (var mstream = new MemoryStream(val))
			using (var reader = new BinaryReader(mstream)){
				int len = reader.ReadInt32();
				var ret = new int[len];
				for (int i=0; i<len; ++i){
					ret[i] = reader.ReadInt32();
				}
				return ret;
			}
		}
		//callback
		public Hybmesh.DCallback callback;
	
		//constructor/destructor
		public Worker(){
			connection = require_connection(hybmesh_exec_path);
			callback = (string s1, string s2, double p1, double p2) => 0;
		}

	};
	public readonly Worker worker;

//public:
	/// <summary> Path to executable (including program file). </summary>
	public static string hybmesh_exec_path =
		"/home/ek/Hybmesh/build/bin/hybmesh";  //>>$EXEPATH
	/// <summary>
	/// Path to shared library libcore_hmconnection_cs (excluding lib file).
	/// Could be changed before initialization of any class instances.
	/// </summary>
	public static string hybmesh_lib_path =
		"/home/ek/Hybmesh/build/lib/core_hybmesh_connection_cs.so";  //>>$LIBPATH

	[Serializable]
	public class EUserInterrupt: Exception{
		public EUserInterrupt(): base ("Interrupted by user"){}
	}
	[Serializable]
	public class ERuntimeError: Exception{
		public ERuntimeError(string message):
			base(message){}
	};
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

	/// <summary> Establishes hybmesh connection </summary>
	public Hybmesh(){
		worker = new Worker();
	}
	protected virtual void Dispose(bool disposing){
		worker.free();
	}
	~Hybmesh(){
		Dispose(false);
	}
	public void Dispose(){
		Dispose(true);
		GC.SuppressFinalize(this);
	}

	//==================== callback
	/// <summary> Callback function delegate.
	/// Returns 0 to continue; 1 to cancel </summary>
	public delegate int DCallback(string s1, string s2, double p1, double p2);

	/// <summary> Assigns non-default callback </summary>
	public void AssignCallback(DCallback cb){
		worker.callback = cb;
	}

	/// <summary> Sets default silent callback </summary>
	public void ResetCallback(){
		worker.callback = (string s1, string s2, double p1, double p2) => 0;
	}

	// =================== Geometrical objects
	public class Object{
		public readonly string sid;
		public readonly Worker worker;
		protected Object(string sid, Worker worker){
			this.sid=sid;
			this.worker=worker;
		}
		//>>$ObjectA
	};
	public class Object2D: Object{
		protected Object2D(string sid, Worker worker): base(sid, worker){}
		//>>$Object2D
	};
	public class Object3D: Object{
		protected Object3D(string sid, Worker worker): base(sid, worker){}
		//>>$Object3D
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
