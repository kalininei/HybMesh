import commongen


class Generator(commongen.Generator):
    def __init__(self):
        super(Generator, self).__init__()

    # =================== virtual methods (could be overwritten)
    @classmethod
    def dict_simple(cls):
        return {
            'TRUE': "true",
            'FALSE': "false",
            'BOOL': "boolean",
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
            'VECBOOL': "boolean[]",
            'VEC_INT_DOUBLE': "Map<Integer, Double>",
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
            'ZEROVECBOOL': 'boolean[]{}',
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
        return 'int[] ' + cls._translate_vec(*args)

    @classmethod
    def _translate_VALVECDOUBLE(cls, *args):
        return 'double[] ' + cls._translate_vec(*args)

    # ====================== abstract methods (should be overwritten)
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

        funcname = cls.to_lower_camel_case(func.name)
        capstring = ["public %s %s(" % (retval, funcname)]

        defargs = []
        func.default_arguments = []
        for a in args:
            if '=' in a:
                a1, a2 = a.split('=')
                a = a1
                if a2 not in ['true', 'false', 'null'] and\
                        a2[0] not in '-.0123456789':
                    defargs.append([a1.split()[1], a2])
                    func.default_arguments.append(['null'] + defargs[-1])
                else:
                    func.default_arguments.append(['val', a1.split()[1], a2])
            capstring.append(a)
            capstring.append(', ')
        if capstring[-1] == ', ':
            capstring.pop()
        capstring.append(')')
        capstring.append(
                ' throws Hybmesh.EUserInterrupt, Hybmesh.ERuntimeError')
        ret = [''.join(capstring)]

        for (k, v) in defargs:
            new = '' if v.startswith('"') else ' new'
            ret.append("{arg} = ({arg} != null) ? {arg} :{new} {val};".format(
                    arg=k, val=v, new=new))
        return ret

    @classmethod
    def _vec_size(cls, vec):
        return "{}.length".format(vec)

    @classmethod
    def _string_init(cls, var, src):
        return "{0} {1} = new {0}({2})".format("String", var, src)

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
        return "if ({0}.length() >= {1}) {0} = {0}.substring({1})".format(
                var, num)

    @classmethod
    def _for_loop(cls, nstring):
        return "for (int i=0; i<{}; ++i)".format(nstring)

    @classmethod
    def _paste_docstring(cls, funccode, func):
        doclines = func.summarystring.split('\n') + ['$'] +\
                func.docstring.split('\n')
        indent = cls.__get_indent(funccode[0])
        dd = [indent + '/**$']
        for line in doclines:
            dd.append(indent + ' * '+line)
        if len(func.default_arguments) > 0:
            dd.append(indent + ' *$')
        for [tp, k, v] in func.default_arguments:
            if tp == 'null':
                dd.append(indent + ' * @param {}'
                          ' if null then {}$'.format(k, v))
            else:
                dd.append(indent + ' * @param {}'
                          ' default is {}$'.format(k, v))

        dd.append(indent + ' */$')
        for i in range(len(dd)):
            funccode.insert(i, dd[i])
