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
            'STRING': "std::string",
        }

    @classmethod
    def dict_simplevec(cls):
        return {
            'VECINT': "std::vector<int>",
            'VECDOUBLE': "std::vector<double>",
            'VECSTRING': "std::vector<std::string>",
            'VECBOOL': "std::vector<bool>",
            'VEC_INT_DOUBLE': "std::vector<std::pair<int, double>>",
        }

    @classmethod
    def dict_hmvec(cls):
        return {k: "std::vector<{}>".format(cls.dict_hm()[k[3:]])
                for k in super(Generator, cls).dict_hmvec().keys()}

    @classmethod
    def dict_zerovec(cls):
        return {k: "{}"
                for k in super(Generator, cls).dict_zerovec().keys()}

    @classmethod
    def dict_none(cls):
        return {
            'NONEPOINT': 'Point2::None()',
            'NONECONTOUR2D': 'Contour2D::None()',
        }

    @classmethod
    def _translate_VALPOINT(cls, *args):
        return cls._translate_VALVECDOUBLE(*args)

    @classmethod
    def _translate_VALPOINT3(cls, *args):
        return cls._translate_VALVECDOUBLE(*args)

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
        ret = 'worker->{}('.format(func)
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

        capstring = ["%s %s(" % (retval, func.name)]
        for a in args:
            capstring.append(a)
            if a is not args[-1]:
                capstring.append(', ')
        capstring.append(')')
        return [''.join(capstring)]

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
    def _vecbyte_init(cls, var, val):
        return "{} {} = {}".format('VecByte', var, val)

    @classmethod
    def _string_pop_begin(cls, var, num):
        return "if ({0}.size() >= {1}) {0} = {0}.substr({1})".format(var, num)

    @classmethod
    def _for_loop(cls, nstring):
        return "for (size_t i=0; i<{}; ++i)".format(nstring)

    @classmethod
    def _paste_docstring(cls, funccode, func):
        doclines = func.summarystring.split('\n') +\
                func.docstring.split('\n')
        indent = cls.__get_indent(funccode[0])
        for i in range(len(doclines)):
            if not doclines[i].endswith('$'):
                doclines[i] = doclines[i] + '$'
            funccode.insert(i, indent + '// '+doclines[i])
