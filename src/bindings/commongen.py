import os


class Generator(object):
    def __init__(self):
        super(Generator, self).__init__()
        # dictionary is filled in subclass constructors
        self.__class__.dictionary = {}
        self.__class__.dictionary.update(self.dict_simple())
        self.__class__.dictionary.update(self.dict_hm())
        self.__class__.dictionary.update(self.dict_simplevec())
        self.__class__.dictionary.update(self.dict_hmvec())
        self.__class__.dictionary.update(self.dict_zerovec())
        self.__class__.dictionary.update(self.dict_none())

        # filled in parse procedure {"OBJECT", "OBJECT2D", ....}
        # containts code that should be pasted to '>>$OBJECT', ...
        self.outstrings = {'Hybmesh': [], 'Contour2D': [], 'Grid2D': [],
                           'Surface3D': [], 'Grid3D': [],
                           'Object2D': [], 'Object3D': [], 'ObjectA': []}

    # ===========  private methods
    @staticmethod
    def __readlines(fn):
        with open(fn) as f:
            return f.readlines()

    @staticmethod
    def __findline(lines, substr):
        for i in range(len(lines)):
            if substr in lines[i]:
                return i
        raise Exception("Substring %s was not found" % substr)

    @staticmethod
    def __get_indent(line):
        s2 = line.strip()
        i = line.find(s2[0])
        return line[:i]

    @staticmethod
    def __substitute(lines, index, sublines):
        del lines[index]
        for line in reversed(sublines):
            lines.insert(index, line)

    # ===================== self methods
    @classmethod
    def _translate(cls, token):
        fnd1 = token.find('(')
        if fnd1 < 0:
            stem = token[1:]
            args = []
        else:
            stem = token[1:fnd1]
            args = token[fnd1 + 1:-1].split(',')

        try:
            return cls.dictionary[stem]
        except:
            try:
                method = getattr(cls, '_translate_%s' % stem)
                return method(*args)
            except:
                raise Exception("undefined stem %s" % stem)

    @classmethod
    def _tos_method(cls, tp, argument):
        if tp == '#STRING':
            fun = '_tos_string'
        elif tp == '#VECPOINT':
            fun = '_tos_vecpoint'
        elif tp == '#VECDOUBLE':
            fun = '_tos_vecdouble'
        elif tp == '#POINT':
            fun = '_tos_point'
        elif tp == '#POINT3':
            fun = '_tos_point3'
        elif tp in ['#GRID2D', '#OBJECT2D', '#CONTOUR2D',
                    '#GRID3D', '#SURFACE3D']:
            fun = '_tos_object'
        elif tp in ['#VECOBJECT2D', '#VECCONTOUR2D', '#VECOBJECT3D',
                    '#VECGRID3D', '#VECGRID2D', '#VECSURFACE3D',
                    '#VECOBJECT']:
            fun = '_tos_vecobject'
        elif tp == '#INT':
            fun = '_tos_int'
        elif tp == '#DOUBLE':
            fun = '_tos_double'
        elif tp == '#BOOL':
            fun = '_tos_bool'
        elif tp == '#VECINT':
            fun = '_tos_vecint'
        elif tp == '#VECSTRING':
            fun = '_tos_vecstring'
        else:
            raise Exception("unknown tos type " + tp)
        return cls._worker_call(fun, argument)

    @classmethod
    def _return_method(cls, tp, arg):
        if tp == "#CONTOUR2D":
            fun = '_to_cont'
        elif tp == "#VECCONTOUR2D":
            fun = '_to_veccont'
        elif tp == "#GRID2D":
            fun = '_to_grid'
        elif tp == "#SURFACE3D":
            fun = '_to_surface'
        elif tp == "#VECSURFACE3D":
            fun = '_to_vecsurface'
        elif tp == "#VECGRID3D":
            fun = '_to_vecgrid3'
        elif tp == "#GRID3D":
            fun = '_to_grid3'
        elif tp == "#VECGRID2D":
            fun = '_to_vecgrid'
        elif tp == "#POINT":
            fun = '_to_point'
        elif tp == "#VECINT":
            fun = '_to_vecint'
        elif tp == '#INT':
            fun = '_to_int'
        elif tp == '#DOUBLE':
            fun = '_to_double'
        else:
            raise Exception("unknown return_method type: " + tp)
        return cls._return_statement(cls._worker_call(fun, arg))

    @classmethod
    def _rawreturn_method(cls, tp, arg):
        if tp == "#VECDOUBLE":
            fun = '_to_vecdouble_raw'
        elif tp == "#VEC_INT_DOUBLE":
            fun = '_to_vec_int_double_raw'
        elif tp == "#VECINT":
            fun = '_to_vecint_raw'
        else:
            raise Exception("unknown rawreturn_method type: " + tp)
        return cls._return_statement(cls._worker_call(fun, arg))

    @classmethod
    def _parse_caption(cls, func):
        define_args = []
        for a in func.args:
            if a[0] == '$ARGHIDDEN':
                continue
            # parse value
            aval = a[-1].split('=')
            s = aval[0]
            assert(len(aval) <= 2)
            if len(aval) == 2:
                if aval[1][0] == '#':
                    s = "%s=%s" % (s, cls._translate(aval[1]))
                else:
                    s = "%s=%s" % (s, aval[1])
            s = "{} {}".format(cls._translate(a[-2]), s)
            s = s.strip()
            define_args.append(s)

        ret = cls._function_caption(define_args, func)
        ret[0] = ret[0] + cls._open_tag()
        return ret

    @classmethod
    def _format_out_string(cls, lines, indent):
        ret = []
        for lline in lines:
            for line in lline:
                if len(line.strip()) > 0 and ';:{}(,.'.find(line[-1]) < 0:
                    line = line + cls._eol_symbol()
                if (len(line) > 0):
                    line = indent + line
                ret.append(line)
        return ret

    @classmethod
    def _parse_args(cls, func):
        subs = []
        args = func.args
        while len(args) > 0:
            a = args.pop(0)
            var = a[-1].split('=')[0]
            if a[0] == '$ARG':
                subs.append(cls._tos_method(a[1], var))
            elif a[0] == '$ARGHIDDEN':
                if a[1][0] == '#':
                    a[1] = cls._translate(a[1])
                subs.append(a[1])
            elif a[0] == '$ARGKW':
                subs.append(cls._concat_strings(
                    '"{}="'.format(a[1]), cls._tos_method(a[2], var)))
            elif a[0].startswith('$ARGSPLIT('):
                alla = [a]
                index = int(a[0][10:-1])
                for c in args:
                    if c[0] == '$ARGSPLIT({})'.format(index):
                        alla.append(c)
                for a in alla[1:]:
                    args.remove(a)
                subs2 = []
                for a2 in alla:
                    var2 = a2[-1].split('=')[0]
                    subs2.append(cls._tos_method(a2[1], var2))
                subs2 = cls._merge_string_code(
                    "tmp{}".format(index), "[]", ", ", subs2, '---')
                subs.extend(subs2)
                subs.append("tmp{}".format(index))
            elif a[0].startswith('$ARGSPLITVEC('):
                alla = [a]
                index = int(a[0][13:-1])
                for c in args:
                    if c[0] == '$ARGSPLITVEC({})'.format(index):
                        alla.append(c)
                for a in alla[1:]:
                    args.remove(a)
                subs2 = ['""',
                         '---' +
                         cls._for_loop(
                             cls._vec_size(alla[0][2].split('=')[0])) +
                         cls._open_tag()]
                for a2 in alla:
                    code = "#" + a2[1][4:]
                    subs2.append(cls._indent() + cls._tos_method(
                         code, a2[2].split('=')[0] + cls._ith(code)))
                if cls._close_tag():
                    subs2.append('---' + cls._close_tag())
                subs2.append('---' + cls._string_pop_begin(
                        'tmp{}'.format(index), 2))
                subs2 = cls._merge_string_code(
                        "tmp{}".format(index), "[]", ", ", subs2, '---')
                subs.extend(subs2)
                subs.append("tmp{}".format(index))
            else:
                raise Exception("unknown argument: " + str(a))
        return subs

    @staticmethod
    def to_upper_camel_case(nm):
        cc = nm.split('_')
        ret = []
        for s in cc:
            lst = list(s)
            lst[0] = lst[0].upper()
            ret.extend(lst)
        ret = ''.join(ret)
        ret = ret.replace('2d', '2D')
        ret = ret.replace('3d', '3D')
        return ret

    @staticmethod
    def to_lower_camel_case(nm):
        ret = Generator.to_upper_camel_case(nm)
        return ret[0].lower() + ret[1:]

    @classmethod
    def _parse_fun(cls, func):
        ret = []
        # parsing arguments
        ret.extend(cls._parse_caption(func))

        # assemble command
        subs = cls._parse_args(func)
        ret.extend(cls._merge_string_code("comstr", "", ", ", subs, ''))

        # return statement
        if len(subs) > 0:
            wc = cls._worker_call(
                    '_apply_command', '"{}"'.format(func.targetfun), 'comstr')
        else:
            wc = cls._worker_call(
                    '_apply_command', '"{}"'.format(func.targetfun), '""')
        if func.argreturn[0] == '$RETURNNO':
            ret.append(wc)
        elif func.argreturn[0] == '$RETURN':
            ret.append(cls._vecbyte_init('comret', wc))
            ret.append(cls._string_init(
                    'sret', cls._worker_call('_tos_vecbyte', 'comret')))
            ret.append(cls._return_method(func.argreturn[1], "sret"))
        elif func.argreturn[0] == '$RETURNRAW':
            ret.append(cls._vecbyte_init('comret', wc))
            ret.append(cls._rawreturn_method(func.argreturn[1], "comret"))
        else:
            raise Exception("invalid return code: %s" % func.argreturn[0])

        # align
        for i in range(1, len(ret)):
            ret[i] = '{}{}'.format(cls._indent(), ret[i])

        # close
        if cls._close_tag():
            ret.append(cls._close_tag())
        ret.append('')
        return ret

    @classmethod
    def _merge_string_code(cls, tmpvar, parant, divider,
                           substrings, starter=''):
        """
        {starter}tmpvar = []
        {starter}tmpvar.append(substrings[0])
        {starter}tmpvar.append(substrings[1])
        {starter}outvar = '[' + {divider}.join(tmpvar) + ']'

        if substring starts with '---',
            simply pastes this substing without '---'.
        if substring has white spaces ahead,
            uses them as string indentation
        """
        ret = []
        string_created = False
        for s in substrings:
            if s.startswith('---'):
                ret.append(s[3:])
            else:
                s2 = s.strip()
                fnd = s.find(s2[0])
                if fnd != 0:
                    indent = s[:fnd]
                else:
                    indent = ''
                if not string_created:
                    ret.append(indent + cls._string_init(tmpvar, s2))
                    string_created = True
                else:
                    ret.append(cls._string_append(indent, tmpvar, s2, ', '))
        if len(parant) == 2:
            ret.append(cls._string_into_parant(tmpvar, parant))

        if starter:
            for i in range(len(ret)):
                ret[i] = "{}{}".format(starter, ret[i])
        return ret

    # =================== virtual methods (could be overwritten)
    @classmethod
    def dict_simple(cls):
        return {
            'TRUE': "",
            'FALSE': "",
            'BOOL': "",
            'INT': "",
            'DOUBLE': "",
            'STRING': "",
        }

    @classmethod
    def dict_hm(cls):
        return {
            'POINT': "Point2",
            'POINT3': "Point3",
            'OBJECT': "Object",
            'OBJECT2D': "Object2D",
            'OBJECT3D': "Object3D",
            'GRID2D': "Grid2D",
            'CONTOUR2D': "Contour2D",
            'GRID3D': "Grid3D",
            'SURFACE3D': "Surface3D",
        }

    @classmethod
    def dict_simplevec(cls):
        return {
            'VECINT': "",
            'VECDOUBLE': "",
            'VECSTRING': "",
            'VECBOOL': "",
            'VEC_INT_DOUBLE': "",
        }

    @classmethod
    def dict_hmvec(cls):
        return {
            'VECPOINT': "",
            'VECPOINT3': "",
            'VECOBJECT': "",
            'VECOBJECT2D': "",
            'VECOBJECT3D': "",
            'VECCONTOUR2D': "",
            'VECGRID2D': "",
            'VECGRID3D': "",
            'VECSURFACE3D': "",
        }

    @classmethod
    def dict_zerovec(cls):
        return {
            'ZEROSTRING': '',
            'ZEROVECOBJECT2D': '',
            'ZEROVECDOUBLE': '',
            'ZEROVECSTRING': '',
            'ZEROVECINT': '',
            'ZEROVECBOOL': '',
            'ZEROVECPOINT': '',
            'ZEROVECPOINT3': '',
            'ZEROVECCONTOUR2D': '',
        }

    @classmethod
    def dict_none(cls):
        return {
            'NONEPOINT': '',
            'NONECONTOUR2D': '',
        }

    @classmethod
    def _write(cls, lines, outfn):
        for i in range(len(lines)):
            if not lines[i].endswith('\n'):
                lines[i] = lines[i] + '\n'
        with open(outfn, 'w') as f:
            f.writelines(lines)

    @classmethod
    def _concat_strings(cls, s1, s2):
        return '{} + {}'.format(s1, s2)

    @classmethod
    def _ith(cls, code):
        return "[i]"

    @classmethod
    def _translate_VALVECINT(cls, *args):
        return cls._translate_vec(*args)

    @classmethod
    def _translate_VALVECDOUBLE(cls, *args):
        return cls._translate_vec(*args)

    @classmethod
    def _translate_VALSTRING(cls, arg):
        return '"%s"' % arg

    @classmethod
    def _translate_VALPOINT(cls, *args):
        return "{}({}, {})".format(cls.dictionary['POINT'], *args)

    @classmethod
    def _translate_VALPOINT3(cls, *args):
        return "{}({}, {}, {})".format(cls.dictionary['POINT3'], *args)

    @classmethod
    def _translate_SID(cls):
        return cls._worker_call('_tos_string', 'sid')

    @classmethod
    def _return_statement(cls, val):
        return "return {}".format(val)

    @classmethod
    def _string_into_parant(cls, s, parant):
        return "{0} = '{1}' + {0} + '{2}'".format(s, parant[0], parant[1])

    @classmethod
    def _indent(cls):
        return '\t'

    @classmethod
    def _eol_symbol(cls):
        return ';'

    @classmethod
    def _close_tag(cls):
        return '}'

    @classmethod
    def _open_tag(cls):
        return '{'

    # ====================== abstract methods (should be overwritten)
    @classmethod
    def _translate_vec(cls, *args):
        raise NotImplementedError

    @classmethod
    def _worker_call(clc, func, *arg):
        raise NotImplementedError

    @classmethod
    def _function_caption(cls, args, func):
        raise NotImplementedError

    @classmethod
    def _vec_size(cls, vec):
        raise NotImplementedError

    @classmethod
    def _string_init(cls, var, src):
        raise NotImplementedError

    @classmethod
    def _string_append(cls, indent, var, what, sep):
        raise NotImplementedError

    @classmethod
    def _vecbyte_init(cls, var, val):
        raise NotImplementedError

    @classmethod
    def _string_pop_begin(cls, var, num):
        raise NotImplementedError

    @classmethod
    def _for_loop(cls, nstring):
        raise NotImplementedError

    # ====================== interface functions
    @staticmethod
    def factory(tp):
        if tp == 'py':
            import pygen
            return pygen.Generator()
        elif tp == 'cpp':
            import cppgen
            return cppgen.Generator()
        elif tp == 'cs':
            import csgen
            return csgen.Generator()
        elif tp == 'java':
            import javagen
            return javagen.Generator()
        elif tp == 'm':
            import mgen
            return mgen.Generator()
        else:
            raise Exception('Unknown type')

    def parse(self, mask):
        """ fills out strings """
        for f in mask.funcs:
            tar = self.outstrings[f.class_]
            tar.append(self._parse_fun(f))

        for v in self.outstrings.values():
            while (v[-1][-1] == ''):
                v[-1].pop()

    def write_to_file(self, basefile, outdir, libpath, exepath):
        # reading basefile and substitution
        bb = self.__readlines(basefile)
        for k, v in self.outstrings.iteritems():
            index = self.__findline(bb, ">>$%s" % k)
            indent = self.__get_indent(bb[index])
            s = self._format_out_string(v, indent)
            self.__substitute(bb, index, s)

        # libdir, exedir
        if libpath is None:
            libpath = "lib-not-defined"
        if exepath is None:
            exepath = "exe-not-defined"
        for i in range(len(bb)):
            if ">>$EXEPATH" in bb[i]:
                quotes = '"' if '"' in bb[i] else "'"
                s3 = bb[i].split(quotes)
                send = '' if s3[2][0] == ' ' else s3[2][0]
                bb[i] = s3[0] + quotes + exepath + quotes + send
            elif ">>$LIBPATH" in bb[i]:
                quotes = '"' if '"' in bb[i] else "'"
                s3 = bb[i].split(quotes)
                send = '' if s3[2][0] == ' ' else s3[2][0]
                bb[i] = s3[0] + quotes + libpath + quotes + send

        # writing to file
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        extension = os.path.splitext(basefile)[1]
        fn = 'Hybmesh' + extension
        self._write(bb, os.path.join(outdir, fn))
