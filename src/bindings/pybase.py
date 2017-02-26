import ctypes as ct
import os
import struct
import ast

class _Worker(object):
    def __init__(self):
        super(_Worker, self).__init__()
        self.connection = -1
        self.callback = lambda n1, n2, p1, p2: 1
        self.cport = None
        self.c_char_data = None
        self.c_char_len = 0

    # low level server communication
    def require_connection(self, path):
        path = path + '\0';
        path = path.encode('utf-8')
        return self.cport.require_connection(path)


    def get_signal(self):
        self.cport.get_signal.restype = ct.c_char;
        a = self.cport.get_signal(ct.c_int(self.connection));
        return a.value;


    def get_data(self):
        if self.c_char_data is not None:
            self.cport.free_char_array(self.c_char_data)
        else:
            self.c_char_data = ct.c_void_p(0);
        self.c_char_len = self.cport.get_data(ct.c_int(self.connection),
                                              ct.byref(self.c_char_data))

        
    def send_signal(self, sig):
        self.cport.send_signal(ct.c_int(self.connection), ct.c_char(sig))


    def send_data(self, sz, data):
        self.cport.send_data(ct.c_int(self.connection), ct.c_int(sz), data)


    def break_connection(self):
        self.cport.break_connection(ct.c_int(self.connection))


    # server communication
    def _send_command(self, func, com):
        self.send_data(len(func), func)
        self.send_data(len(com), com)
        self.send_signal('C')


    def _read_buffer(self):
        self.get_data()
        return self.c_char_len, self.c_char_data


    def _wait_for_signal(self):
        return self.get_signal()


    # command interface
    def _apply_command(self, func, com):
        self._send_command(func, com)
        while 1:
            sig = self._wait_for_signal();
            if sig == 'R':
                return self._read_buffer()
            elif sig == 'B':
                self._apply_callback(self._read_buffer())
            elif sig == 'I':
                raise Hybmesh.EUserInterrupt();
            elif sig == 'E':
                raise Hybmesh.ERuntimeError(self._tos_vecbyte(_read_buffer()));
            elif sig == '0':
                raise Hybmesh.ERuntimeError("Server stopped working");
            else:
                raise Hybmesh.ERuntimeError("Invalid client instruction");

    def _apply_callback(self, buf):
        p1, p2, l1, l2 = struct.unpack_from("=ddii", buf[1])
        n1, n2 = struct.unpack_from("={}s{}s".format(l1, l2), buf[1], 24)
        r = callback(n1, n2, p1, p2)
        sig = 'G' if r == 0 else 'S'
        send_signal(sig)


    # Converters
    def _to_int(self, s):
        return int(s)


    def _to_double(self, s):
        return float(s)


    def _to_point(self, s):
        return ast.literal_eval(s)


    def _to_grid(self, s):
        s = ast.literal_eval(s)
        if s is None:
            return None
        if isinstance(s, list):
            s = s[0]
        return Hybmesh.Grid2D(s, self)


    def _to_cont(self, s):
        s = ast.literal_eval(s)
        if s is None:
            return None
        if isinstance(s, list):
            s = s[0]
        return Hybmesh.Contour2D(s, self)


    def _to_grid3(self, s):
        s = ast.literal_eval(s)
        if s is None:
            return None
        if isinstance(s, list):
            s = s[0]
        return Hybmesh.Grid3D(s, self)


    def _to_surface(self, s):
        s = ast.literal_eval(s)
        if s is None:
            return None
        if isinstance(s, list):
            s = s[0]
        return Hybmesh.Surface3D(s, self)


    def _to_vecint(self, s):
        return ast.literal_eval(s)


    def _to_veccont(self, s):
        s = ast.literal_eval(s)
        if s is None:
            return []
        if not isinstance(s, list):
            s = [s]
        return map(lambda x: Hybmesh.Contour2D(x, self), s)
    

    def _to_vecsurface(self, s):
        s = ast.literal_eval(s)
        if s is None:
            return []
        if not isinstance(s, list):
            s = [s]
        return map(lambda x: Hybmesh.Surface3D(x, self), s)


    def _to_vecgrid(self, s):
        s = ast.literal_eval(s)
        if s is None:
            return []
        if not isinstance(s, list):
            s = [s]
        return map(lambda x: Hybmesh.Grid2D(x, self), s)


    def _to_vecgrid3(self, s):
        s = ast.literal_eval(s)
        if s is None:
            return []
        if not isinstance(s, list):
            s = [s]
        return map(lambda x: Hybmesh.Grid3D(x, self), s)


        #cpp type -> string
        def _tos_bool(self, val):
            return 'True' if val else 'False'

        def _tos_int(self, val):
            return str(val)

        def _tos_double(self, val):
            return "{:.17g}" % val
        
        def _tos_point(self, val):
            if val is not None:
                return "[{:.17g}, {:.17g}]".format(*val)
            else:
                return "None"

        def _tos_point3(self, val):
            if val is not None:
                return "[{:.17g}, {:.17g}, {:.17g}]".format(*val)
            else:
                return "None"

        def _tos_vecbyte(self, val):
            return struct.unpack_from("={}s".format(val[0]), val[1])

        def _tos_vecint(self, val):
            return str(val)

        def _tos_vecdouble(self, val):
            ret = ', '.join(map(lambda x: "{%.17g}".format(x), val))
            return '[' + ret + ']'


        def _tos_vecstring(self, val):
            return str(val)


        def _tos_vecpoint(self, val):
            ret = ', '.join(map(self._tos_point, val))
            return '[' + ret + ']' 


        def _tos_object(self, obj):
            if obj is not None:
                return '"' + obj.sid + '"'
            else:
                return "None"
        
        def _tos_vecobject(self, obj):
            ret = ', '.join(map(self._tos_object), obj)
            return '[' + ret + ']';

        # vecbyte -> cpp type
        def _to_vec_int_double_raw(self, val):
            sz = struct.unpack_from('=i', val[1])
            pos = 4
            ret = []
            for i in range(sz):
                ret.append(list(struct.unpack_from('=id', val, pos)))
                pos += 12
            return ret

        def _to_vecdouble_raw(self, val):
            sz = struct.unpack_from('=i', val[1])
            return list(struct.unpack_from('={}d'.format(sz), val, 4))



class Hybmesh(object):
    hybmesh_exec_path = "home/ek/HybMesh/build/bin/hybmesh"
    hybmesh_lib_path = "home/ek/HybMesh/build/lib/core_hybmesh_connection_py.so"

    class ERuntimeError(Exception):
        def __init__(self, msg):
            super(ERuntimeError, self).__init__("Hybmesh runtime error\n" + msg)

    class EUserInterrupt(Exception):
        def __init__(self):
            super(ERuntimeError, self).__init__("Interrupted by user")

    def __init__(self):
        super(Hybmesh, self).__init__()
        self._worker = _Worker();
        self._worker.connection = self._worker.require_connection(self.hybmesh_exec_path)
        self._worker.cport = ct.cdll.LoadLibrary(os.path.abspath(self.hybmesh_lib_path))

    def __del__(self):
        if self._worker.connection != -1:
            self._worker.break_connection()
            self._worker.connection = -1
        if self._worker.c_char_data is not None:
            self._worker.cport.free_char_array(self._worker.c_char_data)
            self._worker.c_char_data = None

    def assign_callback(self, cb):
            self._worker.callback = cb;

    class Contour2D(object):
        def __init__(self, sid, worker):
            super(Contour2D, self).__init__()
            self.sid = sid
            self._worker = worker
        #>>$Contour2D


    class Grid2D(object):
        def __init__(self, sid, worker):
            super(Grid2D, self).__init__()
            self.sid = sid;
            self._worker = worker;
        #>>$Grid2D


    class Surface3D(object):
        def __init__(self, sid, worker):
            super(Surface3D, self).__init__()
            self.sid = sid;
            self._worker = worker;
        #>>$Surface3D


    class Grid3D(object):
        def __init__(self, sid, worker):
            super(Grid3D, self).__init__()
            self.sid = sid
            self._worker = worker
        #>>$Grid3D

    #>>$Hybmesh
