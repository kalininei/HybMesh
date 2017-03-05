import ctypes as ct
import os
import struct
import ast


class _Worker(object):
    def __init__(self):
        super(_Worker, self).__init__()
        self.callback = lambda n1, n2, p1, p2: 0
        self.cport = ct.cdll.LoadLibrary(os.path.join(
                Hybmesh.hybmesh_lib_path, "libcore_hmconnection_py.so"))
        self.connection = self.require_connection(Hybmesh.hybmesh_exec_path)
        self.c_char_data = None
        self.c_char_len = 0

    def free(self):
        if self.c_char_data is not None:
            self.cport.free_char_array(self._worker.c_char_data)
            self.c_char_data = None
            self.c_char_len = 0
        if self.connection != -1:
            self.break_connection()
            self.connection = -1

    # low level server communication
    def require_connection(self, path):
        path = path + '\0'
        path = path.encode('utf-8')
        return self.cport.require_connection(path)

    def get_signal(self):
        a = self.cport.get_signal(ct.c_int(self.connection))
        return chr(a)

    def get_data(self):
        sz = self.cport.get_data1(ct.c_int(self.connection))
        ret = ct.create_string_buffer(sz)
        self.cport.get_data2(ct.c_int(self.connection), ct.c_int(sz), ret)
        return ret

    def send_signal(self, sig):
        self.cport.send_signal(ct.c_int(self.connection), ct.c_byte(ord(sig)))

    def send_data(self, sz, data):
        data = data.encode('utf-8')
        self.cport.send_data(ct.c_int(self.connection), ct.c_int(sz), data)

    def break_connection(self):
        self.cport.break_connection(ct.c_int(self.connection))

    # server communication
    def _send_command(self, func, com):
        self.send_data(len(func), func)
        self.send_data(len(com), com)
        self.send_signal('C')

    def _read_buffer(self):
        return self.get_data()

    def _wait_for_signal(self):
        return self.get_signal()

    # command interface
    def _apply_command(self, func, com):
        self._send_command(func, com)
        while 1:
            sig = self._wait_for_signal()
            if sig == 'R':
                return self._read_buffer()
            elif sig == 'B':
                self._apply_callback(self._read_buffer())
            elif sig == 'I':
                raise Hybmesh.EUserInterrupt()
            elif sig == 'E':
                raise Hybmesh.ERuntimeError(
                    self._tos_vecbyte(self._read_buffer()))
            elif sig == '0':
                raise Hybmesh.ERuntimeError("Server stopped working")
            else:
                raise Hybmesh.ERuntimeError("Invalid client instruction")

    def _apply_callback(self, buf):
        p1, p2, l1, l2 = struct.unpack_from("=ddii", buf)
        n1, n2 = struct.unpack_from("={}s{}s".format(l1, l2), buf, 24)
        n1, n2 = n1.decode('utf-8'), n2.decode('utf-8')
        r = self.callback(n1, n2, p1, p2)
        sig = 'G' if r == 0 else 'S'
        self.send_signal(sig)

    # Converters
    # removes ending []; splits by ,; strips substrings from """, "'", " ";
    def __parse_vecstring(self, s):
        s = s.strip()
        if (s[0] == '[' and s[-1] == ']'):
            s = s[1:-1]
        if (s == "[]" or s == "None" or len(s) == 0):
            return []
        # split by ',' keeping sublists [x, y] entries untouched
        ret1 = []
        ssymb = '\'\" '
        bracket_level = 0
        for it in map(lambda x: x.strip(ssymb), s.split(',')):
            if len(it) == 0:
                continue
            elif bracket_level == 0:
                ret1.append(it)
            else:
                ret1[-1] = ret1[-1] + ',' + it
            if it[0] == '[':
                bracket_level += 1
            if it[-1] == ']':
                bracket_level -= 1
        return ret1

    def _to_int(self, s):
        return int(s)

    def _to_double(self, s):
        return float(s)

    def _to_point(self, s):
        return ast.literal_eval(s)

    def _to_vecint(self, s):
        return ast.literal_eval(s)

    def _to_object(self, cls, s):
        ss = self.__parse_vecstring(s)
        if (len(ss) == 0 or ss[0] == "None"):
            return None
        else:
            return cls(ss[0], self)

    def _to_vecobject(self, cls, s):
        ss = self.__parse_vecstring(s)
        ret = []
        for sit in ss:
            ret.append(self._to_object(cls, sit))
        return ret

    def _to_grid(self, s):
        return self._to_object(Hybmesh.Grid2D, s)

    def _to_cont(self, s):
        return self._to_object(Hybmesh.Contour2D, s)

    def _to_grid3(self, s):
        return self._to_object(Hybmesh.Grid3D, s)

    def _to_surface(self, s):
        return self._to_object(Hybmesh.Surface3D, s)

    def _to_vecgrid(self, s):
        return self._to_vecobject(Hybmesh.Grid2D, s)

    def _to_veccont(self, s):
        return self._to_vecobject(Hybmesh.Contour2D, s)

    def _to_vecgrid3(self, s):
        return self._to_vecobject(Hybmesh.Grid3D, s)

    def _to_vecsurface(self, s):
        return self._to_vecobject(Hybmesh.Surface3D, s)

    # cpp type -> string
    def _tos_bool(self, val):
        return 'True' if val else 'False'

    def _tos_int(self, val):
        return str(val)

    def _tos_string(self, obj):
        if obj is not None:
            return '"' + obj + '"'
        else:
            return "None"

    def _tos_double(self, val):
        return "{:.16g}".format(val)

    def _tos_point(self, val):
        if val is not None:
            return "[{:.16g}, {:.16g}]".format(*val)
        else:
            return "None"

    def _tos_point3(self, val):
        if val is not None:
            return "[{:.16g}, {:.16g}, {:.16g}]".format(*val)
        else:
            return "None"

    def _tos_vecbyte(self, val):
        return struct.unpack_from(
            "={}s".format(len(val)), val)[0].decode('utf-8')

    def _tos_vecint(self, val):
        return str(val)

    def _tos_vecdouble(self, val):
        if val is None:
            return "None"
        ret = ', '.join(map(lambda x: "{:.16g}".format(x), val))
        return '[' + ret + ']'

    def _tos_vecstring(self, val):
        if val is None:
            return "None"
        return str(val)

    def _tos_vecpoint(self, val):
        if val is None:
            return "None"
        ret = ', '.join(map(self._tos_point, val))
        return '[' + ret + ']'

    def _tos_object(self, obj):
        if obj is not None:
            return self._tos_string(obj.sid)
        else:
            return "None"

    def _tos_vecobject(self, obj):
        if obj is None:
            return "None"
        ret = ', '.join(map(self._tos_object, obj))
        return '[' + ret + ']'

    # vecbyte -> cpp type
    def _to_vec_int_double_raw(self, val):
        sz, = struct.unpack_from('=i', val)
        pos = 4
        ret = []
        for i in range(sz):
            ret.append(list(struct.unpack_from('=id', val, pos)))
            pos += 12
        return ret

    def _to_vecdouble_raw(self, val):
        sz, = struct.unpack_from('=i', val)
        return list(struct.unpack_from('={}d'.format(sz), val, 4))

    def _to_vecint_raw(self, val):
        sz, = struct.unpack_from('=i', val)
        return list(struct.unpack_from('={}i'.format(sz), val, 4))


class Hybmesh(object):
    hybmesh_exec_path = "path/to/hybmesh"  # >>$EXEPATH
    hybmesh_lib_path = "path/to/lib"  # >>$LIBPATH

    class ERuntimeError(Exception):
        def __init__(self, msg):
            super(Hybmesh.ERuntimeError, self).__init__(
                    "Hybmesh runtime error\n" + msg)

    class EUserInterrupt(Exception):
        def __init__(self):
            super(Hybmesh.ERuntimeError, self).__init__(
                    "Interrupted by user")

    def __init__(self):
        super(Hybmesh, self).__init__()
        self._worker = _Worker()

    def __del__(self):
        self._worker.free()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self._worker.free()

    def assign_callback(self, cb):
        self._worker.callback = cb

    def reset_callback(self):
        self._worker.callback = lambda s1, s2, p1, p2: 0

    class Object(object):
        def __init__(self, sid, worker):
            super(Hybmesh.Object, self).__init__()
            self.sid = sid
            self._worker = worker

        # >>$ObjectA

    class Object2D(Object):
        def __init__(self, sid, worker):
            super(Hybmesh.Object2D, self).__init__(sid, worker)

        # >>$Object2D

    class Object3D(Object):
        def __init__(self, sid, worker):
            super(Hybmesh.Object3D, self).__init__(sid, worker)

        # >>$Object3D

    class Contour2D(Object2D):
        def __init__(self, sid, worker):
            super(Hybmesh.Contour2D, self).__init__(sid, worker)

        # >>$Contour2D

    class Grid2D(Object2D):
        def __init__(self, sid, worker):
            super(Hybmesh.Grid2D, self).__init__(sid, worker)

        # >>$Grid2D

    class Surface3D(Object3D):
        def __init__(self, sid, worker):
            super(Hybmesh.Surface3D, self).__init__(sid, worker)

        # >>$Surface3D

    class Grid3D(Object3D):
        def __init__(self, sid, worker):
            super(Hybmesh.Grid3D, self).__init__(sid, worker)

        # >>$Grid3D

    # >>$Hybmesh
