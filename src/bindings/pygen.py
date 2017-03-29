import commongen


class Generator(commongen.Generator):
    def __init__(self):
        super(Generator, self).__init__()

    # =================== virtual methods (could be overwritten)
    @classmethod
    def dict_simple(cls):
        return {
            'TRUE': "True",
            'FALSE': "False",
            'BOOL': "",
            'INT': "",
            'DOUBLE': "",
            'STRING': "",
        }

    @classmethod
    def dict_hm(cls):
        return {k: "" for k in super(Generator, cls).dict_hm().keys()}

    @classmethod
    def dict_zerovec(cls):
        ret = {k: "[]" for k in super(Generator, cls).dict_zerovec().keys()}
        ret["ZEROSTRING"] = '""'
        return ret

    @classmethod
    def dict_none(cls):
        return {
            'NONEPOINT': 'None',
            'NONECONTOUR2D': 'None',
        }

    @classmethod
    def _translate_VALPOINT(cls, *args):
        return cls._translate_VALVECDOUBLE(*args)

    @classmethod
    def _translate_VALPOINT3(cls, *args):
        return cls._translate_VALVECDOUBLE(*args)

    @classmethod
    def _translate_SID(cls):
        return cls._worker_call('_tos_string', 'self.sid')

    @classmethod
    def _indent(cls):
        return '    '

    @classmethod
    def _eol_symbol(cls):
        return ''

    @classmethod
    def _close_tag(cls):
        return ''

    @classmethod
    def _open_tag(cls):
        return ':'

    @classmethod
    def _inline_comment(cls):
        return '#'

    # ====================== abstract methods (should be overwritten)
    @classmethod
    def _translate_vec(cls, *args):
        ret = ['[']
        for a in args:
            ret.append('%s, ' % str(a))
        if len(args) > 0:
            ret[-1] = ret[-1][:-2]
        ret.append(']')
        return ''.join(ret)

    @classmethod
    def _worker_call(clc, func, *arg):
        ret = 'self._worker.{}('.format(func)
        if len(arg) > 0:
            ret = ret + arg[0]
            for i in range(len(arg) - 1):
                ret = ret + ', ' + arg[i+1]
        return ret + ')'

    @classmethod
    def _function_caption(cls, args, func):
        capstring = ["def %s(self, " % func.name]
        for a in args:
            capstring.append(a)
            capstring.append(', ')
        capstring[-1] = capstring[-1][:-2]
        capstring.append(')')
        return [''.join(capstring)]

    @classmethod
    def _vec_size(cls, vec):
        return "len({})".format(vec)

    @classmethod
    def _string_init(cls, var, src):
        return "{} = {}".format(var, src)

    @classmethod
    def _string_append(cls, indent, var, what, sep):
        if sep:
            return "{0}{1} = {1} + '{2}' + {3}".format(indent, var, sep, what)
        else:
            return "{0}{1} = {1} + {2}".format(indent, var, what)

    @classmethod
    def _vecbyte_init(cls, var, val):
        return "{} = {}".format(var, val)

    @classmethod
    def _string_pop_begin(cls, var, num):
        return "{0} = {0}[{1}:]".format(var, num)

    @classmethod
    def _for_loop(cls, nstring):
        return "for i in range({})".format(nstring)

    @classmethod
    def _paste_docstring(cls, funccode, func):
        doclines = []
        indent = cls.__get_indent(funccode[0])
        # summary
        doclines.extend(func.summarystring.split('\n'))
        # arguments types
        for tp, nm in zip(func._argtypes, func._argnames):
            doclines.append("$")
            doclines.append(":param {}: {}".format(nm, cls._doctype(tp)))
        if func.argreturn[0] != "$RETURNNO":
            doclines.append("$")
            doclines.append(":returns: {}".format(
                cls._doctype(func.argreturn[1])))
        # other info
        doclines.append("$")
        doclines.extend(func.docstring.split('\n'))

        doclines.append('"""')
        doclines[0] = '""" ' + doclines[0]
        for i in range(len(doclines)):
            # paste links to hybmeshpack procedure for sphinx documentation
            if "hybmeshpack.hmscript." in doclines[i]:
                istart = doclines[i].find("hybmeshpack.hmscript.")
                iend = doclines[i].find(")", istart)
                bb = doclines[i][istart: iend-1]
                doclines[i] = doclines[i][:istart] + \
                    ":func:`" + bb + '`' + doclines[i][iend+1:]
            funccode.insert(i+1, indent + '    ' + doclines[i])

    # ======================= self methods
    @classmethod
    def _doctype(cls, tp):
        s = ''
        tp = tp[1:]
        if tp.startswith('VEC'):
            s = "list of "
            tp = tp[3:]
        if tp == "_INT_DOUBLE":
            tp = "int, double pairs"
        try:
            tp = {
                "INT": 'int',
                "DOUBLE": 'double',
                "BOOL": 'boolean',
                "STRING": 'string',
                "POINT": '2d point as [x, y]',
                "POINT3": '3d point as [x, y, z]'}[tp]
        except:
            pass
        return s + tp
