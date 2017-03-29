import commongen


class Generator(commongen.Generator):
    def __init__(self):
        super(Generator, self).__init__()

    # ============= virtual overrides
    @classmethod
    def dict_simple(cls):
        return {
            'TRUE': "true",
            'FALSE': "false",
            'BOOL': "bool",
            'INT': "int",
            'DOUBLE': "double",
            'STRING': "String",
        }

    @classmethod
    def dict_simplevec(cls):
        return {
            'VECINT': "int[]",
            'VECDOUBLE': "double[]",
            'VECSTRING': "String[]",
            'VECBOOL': "bool[]",
            'VEC_INT_DOUBLE': "KeyValuePair<int, double>[]",
        }

    @classmethod
    def dict_hmvec(cls):
        return {k: "{}[]".format(cls.dict_hm()[k[3:]])
                for k in super(Generator, cls).dict_hmvec().keys()}

    @classmethod
    def dict_zerovec(cls):
        return {
            'ZEROSTRING': '""',
            'ZEROVECOBJECT2D': 'Object2D[]{}',
            'ZEROVECDOUBLE': 'double[]{}',
            'ZEROVECSTRING': 'String[]{}',
            'ZEROVECINT': 'int[]{}',
            'ZEROVECBOOL': 'bool[]{}',
            'ZEROVECPOINT': 'Point2[]{}',
            'ZEROVECPOINT3': 'Point3[]{}',
            'ZEROVECCONTOUR2D': 'Contour2D[]{}',
        }

    @classmethod
    def dict_none(cls):
        return {
            'NONEPOINT': 'null',
            'NONECONTOUR2D': 'null',
        }

    @classmethod
    def _translate_VALVECINT(cls, *args):
        return 'int[] ' + super(Generator, cls)._translate_VALVECINT(*args)

    @classmethod
    def _translate_VALVECDOUBLE(cls, *args):
        return 'double[] ' + super(Generator, cls)._translate_VALVECINT(*args)

    # ====================== abstract implementations
    @classmethod
    def _translate_vec(cls, *args):
        ret = ['{']
        for a in args:
            ret.append('%s, ' % str(a))
        if len(args) > 0:
            ret[-1] = ret[-1][:-2]
        ret.append('}')
        return ''.join(ret)

    @classmethod
    def _worker_call(clc, func, *arg):
        ret = 'worker.{}('.format(func)
        if len(arg) > 0:
            ret = ret + arg[0]
            for i in range(len(arg) - 1):
                ret = ret + ', ' + arg[i+1]
        return ret + ')'

    @classmethod
    def _function_caption(cls, args, func):
        if func.argreturn[0] == '$RETURNNO':
            retval = 'void'
        else:
            retval = cls._translate(func.argreturn[1])

        funcname = cls.to_upper_camel_case(func.name)
        capstring = ["public %s %s(" % (retval, funcname)]

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
        func.null_argumets = []
        for (k, v) in defargs:
            ret.append("{arg} = {arg} ?? new {val};".format(arg=k, val=v))
            func.null_argumets.append([k, v])
        return ret

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
    def _vecbyte_init(cls, var, val):
        return "{} {} = {}".format('byte[]', var, val)

    @classmethod
    def _string_pop_begin(cls, var, num):
        return "if ({0}.Length >= {1}) {0} = {0}.Substring({1})".format(
                var, num)

    @classmethod
    def _for_loop(cls, nstring):
        return "for (int i=0; i<{}; ++i)".format(nstring)

    @classmethod
    def _paste_docstring(cls, funccode, func):
        doclines = func.summarystring.split('\n') +\
                func.docstring.split('\n')
        indent = cls.__get_indent(funccode[0])
        dd = [indent + '/// <summary>$']
        for line in doclines:
            dd.append(indent + '/// '+line)
        dd.append(indent + '/// </summary>$')
        for [k, v] in func.null_argumets:
            dd.append(indent + '/// <param name="{}">'
                      'if null then {}</param>$'.format(k, v))
        for i in range(len(dd)):
            if not dd[i].endswith('$'):
                dd[i] = dd[i] + '$'
            funccode.insert(i, dd[i])
