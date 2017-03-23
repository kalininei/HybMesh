import commongen
import os


class Generator(commongen.Generator):
    def __init__(self):
        super(Generator, self).__init__()

    @classmethod
    def dict_simple(cls):
        return {
            'TRUE': "true",
            'FALSE': "false",
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
        ret["ZEROSTRING"] = "''"
        return ret

    @classmethod
    def dict_none(cls):
        return {
            'NONEPOINT': 'nan',
            'NONECONTOUR2D': 'nan',
        }

    @classmethod
    def _write(cls, lines, fn):
        # # " -> '
        # for i in range(len(lines)):
        #     lines[i] = lines[i].replace('"', "'")
        # removing superfluous ;, breaking into different files
        for i, ln in enumerate(lines):
            if ln.endswith('\tend;'):
                lines[i] = ln[:-1]
            elif ln.endswith(');') and '\tfunction' in ln:
                lines[i] = ln[:-1]
            elif ln.endswith('];') and '\tfor i' in ln:
                lines[i] = ln[:-1]
        # breaking into Hybmesh.m, Worker.m, ....
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

    @classmethod
    def _concat_strings(cls, s1, s2):
        return "sprintf('%s%s', {}, {})".format(s1, s2)

    @classmethod
    def _ith(cls, code):
        if code in ["#POINT", "#POINT3"]:
            return "(i, :)"
        elif code in ["#CONTOUR2D", "GRID2D", "GRID3D", "SURFACE3D"]:
            return "{i}"
        else:
            return "(i)"

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
    def _translate_SID(cls):
        return cls._worker_call('_tos_string', 'self.sid')

    @classmethod
    def _return_statement(cls, val):
        return "ret = {}".format(val)

    @classmethod
    def _string_into_parant(cls, s, parant):
        return "{0} = sprintf('%s%s%s', '{1}', {0}, '{2}')".format(
                s, parant[0], parant[1])

    @classmethod
    def _close_tag(cls):
        return 'end'

    @classmethod
    def _open_tag(cls):
        return ''

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
        ret = 'self.worker.{}('.format(func)
        if len(arg) > 0:
            ret = ret + arg[0].replace('"', "'")  # fix _apply_command("..")
            for i in range(len(arg) - 1):
                ret = ret + ', ' + arg[i+1].replace('"', "'")
        return ret + ')'

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
    def _vec_size(cls, vec):
        return "size({}, 2)".format(vec)

    @classmethod
    def _string_init(cls, var, src):
        return "{} = {}".format(var, src)

    @classmethod
    def _string_append(cls, indent, var, what, sep):
        if sep:
            return "{0}{1} = sprintf('%s%s%s', {1}, '{2}', {3})".format(
                    indent, var, sep, what)
        else:
            return "{0}{1} = sprintf('%s%s', {1}, {2})".format(
                    indent, var, what)

    @classmethod
    def _vecbyte_init(cls, var, val):
        return "{} = {}".format(var, val)

    @classmethod
    def _string_pop_begin(cls, var, num):
        s = "if length({0})>{1} {0} = {0}({2}: length({0})); end"
        return s.format(var, num, num+1)

    @classmethod
    def _for_loop(cls, nstring):
        return "for i=1:{}".format(nstring)

    @classmethod
    def _paste_docstring(cls, funccode, func):
        doclines = [func.name.upper() + '$']
        indent = cls.__get_indent(funccode[0])
        # summary
        doclines.extend(func.summarystring.split('\n'))
        # arguments types
        for tp, nm in zip(func._argtypes, func._argnames):
            doclines.append("  {}: {}".format(nm, cls._doctype(tp)))
        if func.argreturn[0] != "$RETURNNO":
            doclines.append("  returns: {}".format(
                cls._doctype(func.argreturn[1])))
        # other info
        if func.docstring:
            doclines.append("$")
            doclines.extend(func.docstring.split('\n'))

        funccode.insert(1, '')
        for i in range(len(doclines)):
            funccode.insert(i+1, indent + '% ' + doclines[i])

    # ======================= self methods
    @classmethod
    def _doctype(cls, tp):
        s = ''
        tp = tp[1:]
        if tp == "VEC_INT_DOUBLE":
            tp = "[int, double; ...] pairs"
        if tp.startswith('VEC') and tp[3:8] != "POINT":
            s = "row vector of "
            tp = tp[3:]
        try:
            tp = {
                "INT": 'int',
                "DOUBLE": 'double',
                "BOOL": 'boolean',
                "STRING": 'string',
                "POINT": '2d point as [x, y]',
                "VECPOINT": '2d points list as [x0, y0; x1, y1; ...]',
                "VECPOINT3": '3d points list as [x0, y0, z0; x1, y1, z1; ...]',
                "POINT3": '3d point as [x, y, z]'}[tp]
        except:
            pass
        return s + tp
