import bindparser


class Generator(bindparser.Generator):
    def __init__(self):
        super(Generator, self).__init__()

    dictionary = {
        'TRUE': 'true',
        'FALSE': 'false',
        'ZEROSTRING': '""',
        'ZEROVECOBJECT2D': 'Object2D[]{}',
        'ZEROVECDOUBLE': 'double[]{}',
        'ZEROVECSTRING': 'String[]{}',
        'ZEROVECINT': 'int[]{}',
        'ZEROVECBOOL': 'bool[]{}',
        'ZEROVECPOINT': 'Point2[]{}',
        'ZEROVECPOINT3': 'Point3[]{}',
        'ZEROVECCONTOUR2D': 'Contour2D[]{}',
        'NONEPOINT': 'null',
        'NONECONTOUR2D': 'null',
        'GRID2D': "Grid2D",
        'BOOL': "bool",
        'VECPOINT': "Point2[]",
        'VECPOINT3': "Point3[]",
        'VECINT': "int[]",
        'INT': "int",
        'OBJECT2D': "Object2D",
        'POINT': "Point2",
        'DOUBLE': "double",
        'VECOBJECT2D': "Object2D[]",
        'STRING': "String",
        'VECDOUBLE': "double[]",
        'VECCONTOUR2D': "Contour2D[]",
        'CONTOUR2D': "Contour2D",
        'GRID3D': "Grid3D",
        'VECOBJECT3D': "Object3D",
        'VECGRID2D': "Grid2D[]",
        'VECSURFACE3D': "Surface3D[]",
        'VECGRID3D': "Grid3D[]",
        'VECOBJECT': "Object[]",
        'POINT3': "Point3",
        'VECSTRING': "String[]",
        'VECBOOL': "bool[]",
        'VEC_INT_DOUBLE': "KeyValuePair<int, double>",
    }

    @classmethod
    def __translate_vec(cls, *args):
        ret = ['{']
        for a in args:
            ret.append('%s, ' % str(a))
        if len(args) > 0:
            ret[-1] = ret[-1][:-2]
        ret.append('}')
        return ''.join(ret)

    @classmethod
    def _translate_VALVECINT(cls, *args):
        return "int[]" + cls.__translate_vec(*args)

    @classmethod
    def _translate_VALVECDOUBLE(cls, *args):
        return "double[]" + cls.__translate_vec(*args)

    @classmethod
    def _translate_VALSTRING(cls, arg):
        return '"%s"' % arg

    @classmethod
    def _translate_VALPOINT(cls, *args):
        return "Point2({}, {})".format(*args)

    @classmethod
    def _translate_VALPOINT3(cls, *args):
        return "Point3({}, {}, {})".format(*args)

    @classmethod
    def _translate_SID(cls, arg):
        return cls._worker_call('_tos_string', 'sid');

    @classmethod
    def _worker_call(clc, func, *arg):
        ret = 'worker.{}('.format(func)
        if len(arg) > 0:
            ret = ret + arg[0]
            for i in range(len(arg) - 1):
                ret = ret + ', ' + arg[i+1]
        return ret + ')'

    @classmethod
    def _sid_tos(cls, argument):
        return "{}.sid".format(argument)

    @classmethod
    def _return_statement(cls, val):
        return "return {}".format(val)

    @classmethod
    def _function_caption(cls, args, func):
        if func.argreturn[0] == '$RETURNNO':
            retval = 'void'
        else:
            retval = cls._translate(func.argreturn[1])

        capstring = ["%s %s(" % (retval, func.name)]

        defargs = []
        for a in args:
            if '=' in a:
                a1, a2 = a.split('=')
                if a2 not in ['true', 'false', 'null'] and\
                        a2[0] not in '-.0123456789"':
                    a = a1 + '=null'
                    defargs.append([a1.split()[1], a2])
            capstring.append(a)
            capstring.append(', ')
        if capstring[-1] == ', ':
            capstring.pop()
        capstring.append(')')
        ret = [''.join(capstring)]
        for (k, v) in defargs:
            ret.append("{arg} = {arg} ?? new {val};".format(arg=k, val=v))
        return ret

    @classmethod
    def _close_tag(cls):
        return '}'

    @classmethod
    def _open_tag(cls):
        return '{'

    @classmethod
    def _vec_size(cls, vec):
        return "{}.Length".format(vec)

    @classmethod
    def _string_init(cls, var, src):
        return "String {} = String.Copy({})".format(var, src)

    @classmethod
    def _string_append(cls, indent, var, what, sep):
        if sep:
            return '{0}{1} += "{2}" + {3}'.format(indent, var, sep, what)
        else:
            return '{0}{1} += {2}'.format(indent, var, what)

    @classmethod
    def _string_into_parant(cls, s, parant):
        return "{0} = '{1}' + {0} + '{2}'".format(s, parant[0], parant[1])

    @classmethod
    def _vecbyte_init(cls, var, val):
        return "{} {} = {}".format('byte[]', var, val)

    @classmethod
    def _for_loop(cls, nstring):
        return "for (int i=0; i<{}; ++i)".format(nstring)

    @classmethod
    def _indent(cls):
        return '\t'

    @classmethod
    def _string_pop_begin(cls, var, num):
        return "if ({0}.Length >= {1}) {0} = {0}.Substring({1})".format(
                var, num)

    @classmethod
    def _eol_symbol(cls):
        return ';'
