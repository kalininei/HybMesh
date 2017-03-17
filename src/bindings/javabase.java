import java.lang.reflect.Array;
import java.lang.reflect.Constructor;
import java.util.Arrays;
import java.util.ArrayList;
import java.nio.charset.StandardCharsets;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.Map;
import java.util.HashMap;

class Hybmesh implements AutoCloseable{
	static {
		System.loadLibrary("core_hmconnection_java");
	}
	public static String hybmesh_exec_path = "build/bin/hybmesh";  //>>$EXEPATH

	private static class Worker{
		// ==== communication library
		private int connection = -1;
		private native int require_connection(String path);
		private native byte get_signal(int con);
		private native void send_signal(int con, byte S);
		private native byte[] get_data(int con);
		private native void send_data(int con, byte[] dt);
		private native void break_connection(int con);
	
		//server communication
		private void _send_command(String func, String com){
			byte[] dt1 = func.getBytes(StandardCharsets.UTF_8);
			byte[] dt2 = com.getBytes(StandardCharsets.UTF_8);
			send_signal(connection, (byte)'C');
			send_data(connection, dt1);
			send_data(connection, dt2);
		}
		private byte _wait_for_signal(){
			return get_signal(connection);
		}
		private byte[] _read_buffer(){
			return get_data(connection);
		}
		private void  _apply_callback(byte[] buf){
			ByteBuffer stream = ByteBuffer.wrap(buf);
			stream.order(ByteOrder.nativeOrder());
			double p1 = stream.getDouble();
			double p2 = stream.getDouble();
			int len1 = stream.getInt();
			int len2 = stream.getInt();
			byte[] buf1 = new byte[len1];
			byte[] buf2 = new byte[len2];
			stream.get(buf1, 0, len1);
			stream.get(buf2, 0, len2);
			String n1 = _tos_vecbyte(buf1);
			String n2 = _tos_vecbyte(buf2);
			int r = callback.callback(n1, n2, p1, p2);
			byte R = (r == 0) ? (byte)'G' : (byte)'S';
			send_signal(connection, R);
		}
		public byte[] _apply_command(String func, String com)
				throws Hybmesh.EUserInterrupt, Hybmesh.ERuntimeError{
			_send_command(func, com);
			while (true){
				byte sig = _wait_for_signal();
				switch (sig){
					//normal return
					case 'R': return _read_buffer();
					//callback
					case 'B': _apply_callback(_read_buffer()); break;
					//exception return
					case 'I':
						throw new Hybmesh.EUserInterrupt();
					case 'E':
						throw new Hybmesh.ERuntimeError(
							_tos_vecbyte(_read_buffer()));
					//something went wrong
					case '0':
						throw new Hybmesh.ERuntimeError(
							"Server stopped working");
					default:
						throw new Hybmesh.ERuntimeError(
							"Invalid client instruction");
				}
			}
		}
		//removes ending []; splits by ,; strips substrings from """, "'", " ";
		private String[] __parse_vecstring(String str){
			String s = str.trim();
			if (s.charAt(0)=='[' && s.charAt(s.length()-1)==']')
				s = s.substring(1, s.length()-1);
			if (s == "[]" || s == "None" || s.length() == 0) return new String[]{};
			//split by ',' keeping sublists [x, y] entries untouched
			String[] ret1 = s.split(",");
			for (int i=0; i<ret1.length; ++i){
				int pos1 = 0;
				int pos2 = ret1[i].length();
				while (pos1 < ret1[i].length() &&
				       (ret1[i].charAt(pos1) == ' ' ||
				        ret1[i].charAt(pos1) == '"' ||
				        ret1[i].charAt(pos1) == '\'')) ++pos1;
				while (pos2 > pos1 &&
				       (ret1[i].charAt(pos2-1) == ' ' ||
				        ret1[i].charAt(pos2-1) == '"' ||
				        ret1[i].charAt(pos2-1) == '\'')) --pos2;
				if (pos2 > pos1) ret1[i] = ret1[i].substring(pos1, pos2);
				else ret1[i] = "";
			}
			ArrayList<String> ret2 = new ArrayList<String>();
			int bracket_level=0;
			for (int i=0; i<ret1.length; ++i){
				if (ret1[i].length()==0) continue;
				else if (bracket_level == 0) ret2.add(ret1[i]);
				else ret2.set(
					ret2.size()-1, ret2.get(ret2.size()-1)+','+ret1[i]);
				if (ret1[i].charAt(0) == '[') ++bracket_level;
				if (ret1[i].charAt(ret1[i].length()-1) == ']') --bracket_level;
			}
			String[] ret = new String[ret2.size()];
			return ret2.toArray(ret);
		}
		public int _to_int(String str){
			return Integer.parseInt(str);
		}
		public double _to_double(String str){
			return Double.parseDouble(str);
		}
		public Point2 _to_point(String str){
			if (str == "None") return null;
			String[] s2 = __parse_vecstring(str);
			return new Point2(Double.parseDouble(s2[0]),
					Double.parseDouble(s2[1]));
		}
		private <Obj> Obj _to_object(Class<Obj> type, String str)
				throws Hybmesh.ERuntimeError {
			String[] ss = __parse_vecstring(str);
			if (ss.length == 0 || ss[0] == "None") return null;
			try{
				Constructor<Obj> ctor = type.getDeclaredConstructor(
						String.class, Worker.class);
				return ctor.newInstance(ss[0], this);
			} catch (Exception e){
				throw new Hybmesh.ERuntimeError(e.getMessage());
			}
		}
		private <Obj> Obj[] _to_vecobj(Class<Obj> type, String str)
				throws Hybmesh.ERuntimeError {
			String[] s2 = __parse_vecstring(str);
			@SuppressWarnings("unchecked")
			Obj[] ret = (Obj[]) Array.newInstance(type, s2.length);
			for (int i=0; i<s2.length; ++i){
				ret[i] = _to_object(type, s2[i]);
			}
			return ret;
		}
		public Grid2D _to_grid(String str) throws Hybmesh.ERuntimeError{
			return _to_object(Grid2D.class, str);
		}
		public Contour2D _to_cont(String str) throws Hybmesh.ERuntimeError{
			return _to_object(Contour2D.class, str);
		}
		public Grid3D _to_grid3(String str) throws Hybmesh.ERuntimeError{
			return _to_object(Grid3D.class, str);
		}
		public Surface3D _to_surface(String str) throws Hybmesh.ERuntimeError{
			return _to_object(Surface3D.class, str);
		}
		public int[] _to_vecint(String str){
			String[] s2 = __parse_vecstring(str);
			int[] ret = new int[s2.length];
			for (int i=0; i<ret.length; ++i){
				ret[i] = Integer.parseInt(s2[i]);
			}
			return ret;
		}
		public Contour2D[] _to_veccont(String str) throws Hybmesh.ERuntimeError{
			return _to_vecobj(Contour2D.class, str);
		}
		public Grid2D[] _to_vecgrid(String str) throws Hybmesh.ERuntimeError{
			return _to_vecobj(Grid2D.class, str);
		}
		public Surface3D[] _to_vecsurface(String str) throws Hybmesh.ERuntimeError{
			return _to_vecobj(Surface3D.class, str);
		}
		public Grid3D[] _to_vecgrid3(String str) throws Hybmesh.ERuntimeError{
			return _to_vecobj(Grid3D.class, str);
		}
		// cpp type -> string
		public String _tos_bool(boolean val){
			return val ? "True" : "False";
		}
		public String _tos_int(int val){
			return Integer.toString(val);
		}
		public String _tos_double(double val){
			return Double.toString(val);
		}
		public String _tos_string(String val){
			if (val == null) return "None";
			return '"' + val + '"';
		}
		public String _tos_point(Point2 val){
			if (val == null) return "None";
			return "[" + Double.toString(val.x) + ", " +
				Double.toString(val.y) + "]";
		}
		public String _tos_point3(Point3 val){
			if (val == null) return "None";
			return "[" + Double.toString(val.x) + ", " +
				Double.toString(val.y) + ", " +
				Double.toString(val.z) + "]";
		}
		public String _tos_vecbyte(byte[] val){
			if (val == null) return "None";
			return new String(val, StandardCharsets.UTF_8);
		}
		public String _tos_vecint(int[] val){
			if (val == null) return "None";
			return Arrays.toString(val);
		}
		public String _tos_vecdouble(double[] val){
			if (val == null) return "None";
			return Arrays.toString(val);
		}
		public String _tos_vecstring(String[] val){
			if (val == null) return "None";
			String[] sval = new String[val.length];
			for (int i=0; i<val.length; ++i){
				sval[i] = '"' + val[i] + '"';
			}
			return Arrays.toString(sval);
		}
		public String _tos_vecpoint(Point2[] val){
			if (val == null) return "None";
			String[] s = new String[val.length];
			for (int i=0; i<val.length; ++i){
				s[i] = _tos_point(val[i]);
			}
			return '[' + String.join(", ", s) + ']';
		}
		public String _tos_object(Hybmesh.Object val){
			if (val == null) return "None";
			else return '"' + val.sid + '"';
		}
		public String _tos_vecobject(Hybmesh.Object[] val){
			if (val == null) return "None";
			if (val.length == 0) return "[]";
			String ret = new String();
			for (int i=0; i<val.length; ++i){
				ret += _tos_object(val[i]) + ", ";
			}
			return '[' + ret + ']';
		}
		// vecbyte -> cpp type
		public Map<Integer, Double> _to_vec_int_double_raw(byte[] val){
			ByteBuffer stream = ByteBuffer.wrap(val);
			stream.order(ByteOrder.nativeOrder());
			int sz = stream.getInt();
			Map<Integer, Double> ret = new HashMap<Integer, Double>();
			for (int i=0; i<sz; ++i){
				int a = stream.getInt();
				double b = stream.getDouble();
				ret.put(a, b);
			}
			return ret;
		}
		public int[] _to_vecint_raw(byte[] val){
			ByteBuffer stream = ByteBuffer.wrap(val);
			stream.order(ByteOrder.nativeOrder());
			int sz = stream.getInt();
			int[] ret = new int[sz];
			for (int i=0; i<sz; ++i){
				ret[i] = stream.getInt();
			}
			return ret;
		}
		public double[] _to_vecdouble_raw(byte[] val){
			ByteBuffer stream = ByteBuffer.wrap(val);
			stream.order(ByteOrder.nativeOrder());
			int sz = stream.getInt();
			double[] ret = new double[sz];
			for (int i=0; i<sz; ++i){
				ret[i] = stream.getDouble();
			}
			return ret;
		}
		//callback
		public Hybmesh.ICallback callback;
		public static class DefaultCallback implements Hybmesh.ICallback{
			public int callback(String n1, String n2, double p1, double p2){
				return 0;
			}
		};
		//constructor/destructor
		public Worker(){
			connection = require_connection(hybmesh_exec_path);
			callback = new DefaultCallback();
		}
		public void dispose(){
			if (connection != -1){
				break_connection(connection);
			}
			connection = -1;
		}
		@Override
		protected void finalize() throws Throwable{
			try{
				dispose();
			} finally {
				super.finalize();
			}
		}
	};
	private final Worker worker;

//public:
	public Hybmesh(){
		worker=new Worker();
	};
	@Override
	public void close(){
		worker.dispose();
	}

	public static class ERuntimeError extends Exception{
		public static final long serialVersionUID = 1L;
		public ERuntimeError(String message){
			super(message);
		}
	};
	public static class EUserInterrupt extends Exception{
		public static final long serialVersionUID = 1L;
		public EUserInterrupt(){
			super("Interrupted by user");
		}
	};

	/**
	 * Callback function.
	 *
	 * Should return 0 to proceed with operation, 1 for cancellation request.
	 */
	public interface ICallback{
		public int callback(String n1, String n2, double p1, double p2);
	}

	/**
	 * Assigns non-default callback.
	 */
	public void assignCallback(ICallback cb){
		worker.callback = cb;
	}

	/**
	 * Sets default silent callback.
	 */
	public void resetCallback(){
		worker.callback = new Worker.DefaultCallback();
	}

	public static class Point2{
		public double x, y;
		public Point2(double x, double y){
			this.x = x;
			this.y = y;
		}
	};
	public static class Point3{
		public double x, y, z;
		public Point3(double x, double y, double z){
			this.x = x;
			this.y = y;
			this.z = z;
		}
	};
	public static class Object{
		public final String sid;
		protected final Worker worker;
		public Object(String sid_, Worker worker_){
			sid = sid_;
			worker = worker_;
		}
		//>>$ObjectA
	};
	public static class Object2D extends Object{
		public Object2D(String sid, Worker worker){
			super(sid, worker);
		}
		//>>$Object2D
	};
	public static class Object3D extends Object{
		public Object3D(String sid, Worker worker){
			super(sid, worker);
		}
		//>>$Object3D
	};
	public static class Contour2D extends Object2D{
		public Contour2D(String sid, Worker worker){
			super(sid, worker);
		};
		//>>$Contour2D
	};
	public static class Grid2D extends Object2D{
		public Grid2D(String sid, Worker worker){
			super(sid, worker);
		};
		//>>$Grid2D
	};
	public static class Surface3D extends Object3D{
		public Surface3D(String sid, Worker worker){
			super(sid, worker);
		}
		//>>$Surface3D
	};

	public static class Grid3D extends Object3D{
		public Grid3D(String sid, Worker worker){
			super(sid, worker);
		}
		//>>$Grid3D
	};

	//>>$Hybmesh
};
