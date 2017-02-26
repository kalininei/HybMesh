import bindparser
import os

class Generator(bindparser.Generator):
    def __init__(self):
        super(Generator, self).__init__()

    dictionary = {
        'TRUE': 'true',
        'FALSE': 'false',
        'ZEROSTRING': '""',
        'ZEROVECOBJECT2D': '[]',
        'ZEROVECDOUBLE': '[]',
        'ZEROVECSTRING': '[]',
        'ZEROVECINT': '[]',
        'ZEROVECBOOL': '[]',
        'ZEROVECPOINT': '[]',
        'ZEROVECPOINT3': '[]',
        'ZEROVECCONTOUR2D': '[]',
        'NONEPOINT': 'nan',
        'NONECONTOUR2D': 'nan',
        'GRID2D': "",
        'BOOL': "",
        'VECPOINT': "",
        'VECPOINT3': "",
        'VECINT': "",
        'INT': "",
        'OBJECT2D': "",
        'POINT': "",
        'DOUBLE': "",
        'VECOBJECT2D': "",
        'STRING': "",
        'VECDOUBLE': "",
        'VECCONTOUR2D': "",
        'CONTOUR2D': "",
        'GRID3D': "",
        'VECOBJECT3D': "",
        'VECGRID2D': "",
        'VECSURFACE3D': "",
        'VECGRID3D': "",
        'VECOBJECT': "",
        'POINT3': "",
        'VECSTRING': "",
        'VECBOOL': "",
        'VEC_INT_DOUBLE': "",
    }


    @classmethod
    def __translate_vec(cls, *args):
        ret = ['[']
        for a in args:
            ret.append('%s, ' % str(a))
        if len(args) > 0:
            ret[-1] = ret[-1][:-2]
        ret.append(']')
        return ''.join(ret)


    @classmethod
    def _translate_VALVECINT(cls, *args):
        return cls.__translate_vec(*args)


    @classmethod
    def _translate_VALVECDOUBLE(cls, *args):
        return cls.__translate_vec(*args)


    @classmethod
    def _translate_VALSTRING(cls, arg):
        return "'%s'" % arg


    @classmethod
    def _translate_VALPOINT(cls, *args):
        return cls._translate_VALVECDOUBLE(*args)


    @classmethod
    def _translate_VALPOINT3(cls, *args):
        return cls._translate_VALVECDOUBLE(*args)


    @classmethod
    def _translate_SELF(cls, arg):
        return "self.{}".format(arg)

    @classmethod
    def _worker_call(clc, func, *arg):
        ret = 'self.worker.{}('.format(func)
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
        return "ret = {}".format(val)
    
    @classmethod
    def _function_caption(cls, args, func):
        capstring = ["function ret=%s(self, " % func.name]
        for a in args:
            capstring.append(a)
            capstring.append(', ')
        capstring[-1] = capstring[-1][:-2]
        capstring.append(')')
        return [''.join(capstring)]

    @classmethod
    def _close_tag(cls):
        return 'end'

    @classmethod
    def _open_tag(cls):
        return ''

    @classmethod
    def _vec_size(cls, vec):
        return "size({}, 2)".format(vec)

    @classmethod
    def _string_init(cls, var, src):
        return "{} = {}".format(var, src)

    @classmethod
    def _string_append(cls, indent, var, what, sep):
        if sep:
            return '{0}{1} = {1} + "{2}" + {3}'.format(indent, var, sep, what)
        else:
            return "{0}{1} = {1} + {2}".format(indent, var, what)

    @classmethod
    def _string_into_parant(cls, s, parant):
        return '{0} = "{1}" + {0} + "{2}"'.format(s, parant[0], parant[1])

    @classmethod
    def _vecbyte_init(cls, var, val):
        return "{} = {}".format(var, val)

    @classmethod
    def _for_loop(cls, nstring):
        return "for i in [1:{}]".format(nstring)

    @classmethod
    def _indent(cls):
        return '\t'

    @classmethod
    def _string_pop_begin(cls, var, num):
        return "if size({0}, 2)>{1} {0} = {0}(num + 1: size({0}, 2)) end".format(var, num)

    @classmethod
    def _eol_symbol(cls):
        return ';'

    @classmethod
    def _write(cls, lines, fn):
        "removing superfluous ;, breaking into different files"
        for i, ln in enumerate(lines):
            if ln.endswith('\tend;'):
                lines[i] = ln[:-1]
            elif ln.endswith(');') and '\tfunction' in ln:
                lines[i] = ln[:-1]
            elif ln.endswith('];') and '\tfor i' in ln:
                lines[i] = ln[:-1]
        "breaking into Hybmesh.m, Worker.m, ...."
        istart = []
        for i, ln in enumerate(lines):
            if ln.startswith('classdef'):
                istart.append(i)
        istart.append(len(lines))
        for i in range(len(istart)-1):
            classname = lines[istart[i]].split()[1]
            bb = lines[istart[i]:istart[i+1]]
            mfn = os.path.join(os.path.split(fn)[0], classname + '.m')
            super(Generator, cls)._write(bb, mfn)
