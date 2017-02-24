"""Script for assembling binding headers.
For installation purposes only"""
import sys

class Func(object):
    def __init__(self):
        self.class_ = None
        self.name = None
        self.targetfun = None
        self.args = []
        self.argreturn = None

    def supplement(self):
        if self.name is None:
            raise Exception("function without name")
        if self.argreturn is None:
            raise Exception("function without return value: " + self.name)
        if self.class_ is None:
            self.class_ = 'Hybmesh'
        if self.targetfun is None:
            self.targetfun = self.name


class MaskParser(object):
    def __init__(self, fn):
        super(MaskParser, self).__init__()
        self.funcs = []
        strfun = self.__extract_funcs(fn)
        for sf in strfun:
            f = Func()
            tok = self.__tokenize(sf)
            for it in tok:
                t3 = it[0][1:4]
                if t3 == 'ARG':
                    f.args.append(it)
                elif t3 == 'RET':
                    f.argreturn = it
                elif t3 == 'NAM':
                    f.name = it[1]
                elif t3 == 'CLA':
                    f.class_ = it[1]
                elif t3 == 'TFU':
                    f.targetfun = it[1]
            try:
                f.supplement()
            except:
                print tok
                raise
            self.funcs.append(f)


    @staticmethod
    def __extract_funcs(fn):
        ret = []
        with open(fn) as f:
            stream = f.read()
            icur = 0
            while 1:
                istart = stream.find('$FUNC(', icur)
                if istart < 0:
                    break
                icur = istart + 6
                iend = stream.find('$FUNC', icur)
                if iend < 0 or stream[iend + 5] != ')':
                    print stream[istart:iend]
                    raise Exception('no matching $FUNC)')
                icur = iend
                substr = stream[istart+6:iend]
                ret.append(substr.strip())
        return ret


    @staticmethod
    def __tokenize(lines):
        ret = []
        words = lines.split()
        toks = []
        for i, w in enumerate(words):
            if w.startswith('$'):
                toks.append(i)
        toks.append(len(words))
        for i, itok in enumerate(toks[:-1]):
            ret.append(words[itok: toks[i+1]])
        return ret




class Generator(object):
    def __init__(self):
        super(Generator, self).__init__()
        self.outstrings = {}

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

    # =========== for overriding
    def _format_out_string(self, lines, indent):
        raise NotImplementedError


    def _parse_fun(self, func):
        "-> [str]"
        raise NotImplementedError


    # ===========  self methods
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


    @staticmethod
    def __write(lines, outfn):
        for i in range(len(lines)):
            if not lines[i].endswith('\n'):
                lines[i] = lines[i] + '\n'
        with open(outfn, 'w') as f:
            f.writelines(lines)


    @classmethod
    def _translate(cls, token):
        fnd1 = token.find('(')
        if fnd1 < 0:
            stem = token[1:]
            args = []
        else:
            stem = token[1:fnd1]
            args = token[fnd1 + 1 : -1].split(',')

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
        if tp in ['#GRID2D', '#OBJECT2D', '#CONTOUR2D', '#GRID3D']:
            return cls._sid_tos(argument)
        elif tp == '#STRING':
            return argument
        else:
            if tp == '#VECPOINT':
                fun = '_tos_vecpoint'
            elif tp == '#VECDOUBLE':
                fun = '_tos_vecdouble'
            elif tp == '#POINT':
                fun = '_tos_point'
            elif tp == '#POINT3':
                fun = '_tos_point3'
            elif tp in ['#VECOBJECT2D', '#VECCONTOUR2D', '#VECOBJECT3D', '#VECGRID3D',
                        '#VECGRID2D', '#VECSURFACE3D', '#VECOBJECT']:
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
        elif tp == "#MAP_INT_DOUBLE":
            fun = '_to_map_int_double_raw'
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
                    s = "%s %s=%s" % (cls._translate(a[-2]), s, cls._translate(aval[1]))
                else:
                    s = "%s %s=%s" % (cls._translate(a[-2]), s, aval[1])
            s = s.strip()
            define_args.append(s)

        return cls._function_caption(define_args, func)


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
                subs.append("'{}=' + {}".format(a[1], cls._tos_method(a[2], var)))
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
                subs2 = cls._merge_string_code("tmp{}".format(index), "[]", ", ", subs2, '---')
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
                         '---' + cls._for_loop(cls._vec_size(alla[0][2].split('=')[0]))]
                for a2 in alla:
                    subs2.append(cls._indent() + cls._tos_method(
                        "#" + a2[1][4:], a2[2].split('=')[0] + '[i]'))
                if cls._close_tag():
                    subs2.append('---' + cls._close_tag())
                subs2.append('---' + cls._string_pop_begin('tmp{}'.format(index), 2))
                subs2 = cls._merge_string_code("tmp{}".format(index), "[]", ", ", subs2, '---')
                subs.extend(subs2)
                subs.append("tmp{}".format(index))
            else:
                raise Exception("unknown argument: " + str(a))
        return subs


    @classmethod
    def _parse_fun(cls, func):
        ret = []
        # parsing arguments
        ret.append(cls._parse_caption(func))

        # assemble command
        subs = cls._parse_args(func)
        ret.extend(cls._merge_string_code("comstr", "()", ", ", subs, ''))

        # return statement
        wc = cls._worker_call('_apply_command', func.targetfun, 'comstr')
        if func.argreturn[0] == '$RETURNNO':
            ret.append(wc)
        elif func.argreturn[0] == '$RETURN':
            ret.append(cls._vecbyte_init('comret', wc))
            ret.append(cls._string_init('sret', cls._worker_call('_tos_vecbyte', 'comret')))
            ret.append(cls._return_method(func.argreturn[1], "sret"))
        elif func.argreturn[0] == '$RETURNRAW':
            ret.append(cls._vecbyte_init('comret', wc))
            ret.append(cls._rawreturn_method(func.argreturn[1], "comret"))
        else:
            raise Exception("invalid return code: %s" % func.argreturn[0])

        # align
        for i in range(1, len(ret)):
            ret[i] = '{}{}'.format(cls._indent(), ret[i])
        return ret
    

    @classmethod
    def _merge_string_code(cls, tmpvar, parant, divider, substrings, starter=''):
        """
        {starter}tmpvar = []
        {starter}tmpvar.append(substrings[0])
        {starter}tmpvar.append(substrings[1])
        {starter}outvar = '[' + {divider}.join(tmpvar) + ']'

        if substring starts with '---', simply pastes this substing without '---'.
        if substring has white spaces ahead, uses them as string indentation
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
        ret.append(cls._string_into_parant(tmpvar, parant))

        if starter:
            for i in range(len(ret)):
                ret[i] = "{}{}".format(starter, ret[i])
        return ret





        
    # interface functions
    def parse(self, mask):
        """ fills out strings """
        self.outstrings = {'Hybmesh':[], 'Contour2D': [], 'Grid2D': [],
                           'Surface3D': [], 'Grid3D': []}
        for f in mask.funcs:
            tar = self.outstrings[f.class_]
            tar.append(self._parse_fun(f))


    def write_to_file(self, basefile, outfile):
        bb = self.__readlines(basefile)
        for k, v in self.outstrings.iteritems():
            index = self.__findline(bb, ">>$%s" % k)
            indent = self.__get_indent(bb[index])
            s = self._format_out_string(v, indent)
            self.__substitute(bb, index, s)
        self.__write(bb, outfile)



if __name__ == "__main__":
    """arguments:
        -target py/cpp/cs/m/java
        -mask funcmask_file
        -base basefile with >>$ instructions
        -out output_file
    """
    sys.argv = ['-target', 'py', '-mask', 'funcmask', '-base', 'pybase.py', '-out', 'out.py']
    target, mask, base, out = None, None, None, None
    for i, nm in enumerate(sys.argv):
        if nm == '-target':
            target = sys.argv[i + 1];
        elif nm == '-mask':
            mask = sys.argv[i + 1];
        elif nm == '-base':
            base = sys.argv[i + 1];
        elif nm == '-out':
            out = sys.argv[i + 1];
    if None in [target, mask, base, out]:
        raise Exception("Incorrect arguments")

    mask = MaskParser(mask)
    generator = Generator.factory(target)
    generator.parse(mask)
    generator.write_to_file(base, out)
