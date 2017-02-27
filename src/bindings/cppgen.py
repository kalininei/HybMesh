import bindparser


class Generator(bindparser.Generator):
    def __init__(self):
        super(Generator, self).__init__()

    dictionary = {
        'TRUE': 'true',
        'FALSE': 'false',
        'ZEROSTRING': '""',
        'ZEROVECOBJECT2D': 'std::vector<Object2D>()',
        'ZEROVECDOUBLE': 'std::vector<double>()',
        'ZEROVECSTRING': 'std::vector<std::string>()',
        'ZEROVECINT': 'std::vector<int>()',
        'ZEROVECBOOL': 'std::vector<bool>()',
        'ZEROVECPOINT': 'std::vector<Point2>()',
        'ZEROVECPOINT3': 'std::vector<Point3>()',
        'ZEROVECCONTOUR2D': 'std::vector<Contour2D>()',
        'NONEPOINT': 'Point2::None()',
        'NONECONTOUR2D': 'Contour2D::None()',
        'GRID2D': "Grid2D",
        'BOOL': "bool",
        'VECPOINT': "std::vector<Point2>",
        'VECPOINT3': "std::vector<Point3>",
        'VECINT': "std::vector<int>",
        'INT': "int",
        'OBJECT2D': "Object2D",
        'POINT': "Point2",
        'DOUBLE': "double",
        'VECOBJECT2D': "std::vector<Object2D>",
        'STRING': "std::string",
        'VECDOUBLE': "std::vector<double>",
        'VECCONTOUR2D': "std::vector<Contour2D>",
        'CONTOUR2D': "Contour2D",
        'GRID3D': "Grid3D",
        'VECOBJECT3D': "std::vector<Object3D>",
        'VECGRID2D': "std::vector<Grid2D>",
        'VECSURFACE3D': "std::vector<Surface3D>",
        'VECGRID3D': "std::vector<Grid3D>",
        'VECOBJECT': "std::vector<Object>",
        'POINT3': "Point3",
        'VECSTRING': "std::vector<std::string>",
        'VECBOOL': "std::vector<bool>",
        'VEC_INT_DOUBLE': "std::vector<std::pair<int, double> >",
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
        return cls.__translate_vec(*args)

    @classmethod
    def _translate_VALVECDOUBLE(cls, *args):
        return cls.__translate_vec(*args)

    @classmethod
    def _translate_VALSTRING(cls, arg):
        return '"%s"' % arg

    @classmethod
    def _translate_VALPOINT(cls, *args):
        return cls._translate_VALVECDOUBLE(*args)

    @classmethod
    def _translate_VALPOINT3(cls, *args):
        return cls._translate_VALVECDOUBLE(*args)

    @classmethod
    def _translate_SELF(cls, arg):
        return "{}".format(arg)

    @classmethod
    def _worker_call(clc, func, *arg):
        ret = 'worker->{}('.format(func)
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
        for a in args:
            capstring.append(a)
            if a is not args[-1]:
                capstring.append(', ')
        capstring.append(')')
        return [''.join(capstring)]

    @classmethod
    def _close_tag(cls):
        return '}'

    @classmethod
    def _open_tag(cls):
        return '{'

    @classmethod
    def _vec_size(cls, vec):
        return "{}.size()".format(vec)

    @classmethod
    def _string_init(cls, var, src):
        return "{} {} = {}".format("std::string", var, src)

    @classmethod
    def _string_append(cls, indent, var, what, sep):
        if sep:
            return '{0}{1} += std::string("{2}") + {3}'.format(
                    indent, var, sep, what)
        else:
            return '{0}{1} += {2}'.format(indent, var, what)

    @classmethod
    def _string_into_parant(cls, s, parant):
        return "{0} = '{1}' + {0} + '{2}'".format(s, parant[0], parant[1])

    @classmethod
    def _vecbyte_init(cls, var, val):
        return "{} {} = {}".format('VecByte', var, val)

    @classmethod
    def _for_loop(cls, nstring):
        return "for (int i=0; i<{}; ++i)".format(nstring)

    @classmethod
    def _indent(cls):
        return '\t'

    @classmethod
    def _string_pop_begin(cls, var, num):
        return "if ({0}.size() >= {1}) {0} = {0}.substr({1})".format(var, num)

    @classmethod
    def _eol_symbol(cls):
        return ';'
