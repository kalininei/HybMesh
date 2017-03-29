"""Script for assembling binding headers.
For installation purposes only"""
import sys
import inspect
import commongen


class Func(object):
    def __init__(self):
        self.class_ = None
        self.name = None
        self.targetfun = None
        self.args = []
        self.argreturn = None
        self.summarystring = None
        self.docstring = None
        self._argnames = []
        self._argtypes = []
        self._argdefaults = []

    def supplement(self):
        if self.name is None:
            raise Exception("function without name")
        if self.argreturn is None:
            raise Exception("function without return value: " + self.name)
        if self.class_ is None:
            self.class_ = 'Hybmesh'
        if self.targetfun is None:
            self.targetfun = self.name

        # documentation data
        if self.docstring is None:
            import hybmeshpack.hmscript as hm
            if '_bindoper' in self.targetfun:
                f = getattr(hm._bindoper, self.targetfun[10:])
            else:
                f = getattr(hm, self.targetfun)
            if inspect.getdoc(f) is not None:
                doc = inspect.getdoc(f)
                isent1, isent2 = doc.find('.'), doc.find('\n\n')
                if isent1 < 0:
                    isent1 = isent2
                if isent2 < 0:
                    isent2 = isent1
                isent = min(isent1, isent2)
                self.docstring = doc[:isent] + '.'
            else:
                self.docstring = "not documented yet."
        if '_bindoper' not in self.targetfun:
            ainfo = "See details in " +\
                "hybmeshpack.hmscript.{}().".format(self.targetfun)
            self.docstring = self.docstring + '\n' + ainfo
        dd = self.docstring.strip().split('.', 1)
        self.summarystring = (dd[0]+".").strip()
        if len(dd) > 1:
            self.docstring = dd[1].strip()
            self.docstring = self.docstring.replace('\n', '$\n')
        else:
            self.docstring = '$'
        self.summarystring = self.summarystring.replace('\n', '$\n')

        # _argnames, types, defaults
        for a in self.args:
            if a[0] != "$ARGHIDDEN":
                self._argtypes.append(a[-2])
                sp = a[-1].split('=')
                self._argnames.append(sp[0])
                sp2 = sp[1] if len(sp) > 1 else None
                self._argdefaults.append(sp2)


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
                elif t3 == 'DOC':
                    f.docstring = it[1]
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
        # here we remove $DOC section because it should not
        # be splitted by whitespaces like all others.
        docstart = lines.find("$DOC")
        if docstart >= 0:
            docend = lines.find("$", docstart+4)
            docstr = lines[docstart+4:docend]
            docstr = docstr.replace("\t", "").strip()
            lines = lines[:docstart] + lines[docend:]

        words = lines.split()
        toks = []
        for i, w in enumerate(words):
            if w.startswith('$'):
                toks.append(i)
        toks.append(len(words))
        for i, itok in enumerate(toks[:-1]):
            ret.append(words[itok: toks[i+1]])

        # bring not splitted $DOC section back
        if docstart >= 0:
            ret.append(["$DOC", docstr])

        return ret


if __name__ == "__main__":
    """possible arguments:
        [py/cpp/java/m/oct/all]
        [hmdir /path/to/hybmeshpack/directory]
        [odir /path/to/output/directory]
        [libpath /path/to/lib/]
        [exepath /path/to/exe/]
        [vers '0.1.2']
    """
    try:
        dirindex = sys.argv.index('odir')
        odir = sys.argv[dirindex + 1]
    except ValueError:
        odir = 'out'

    try:
        dirindex = sys.argv.index('hmdir')
        hmdir = sys.argv[dirindex + 1]
        sys.path.insert(0, hmdir)
    except ValueError:
        raise Exception("no hmdir defined")

    mask = 'funcmask'

    if 'libpath' in sys.argv:
        ind = sys.argv.index('libpath')
        libpath = sys.argv[ind+1]
    else:
        libpath = None

    if 'exepath' in sys.argv:
        ind = sys.argv.index('exepath')
        exepath = sys.argv[ind+1]
    else:
        exepath = None

    if 'vers' in sys.argv:
        ind = sys.argv.index('vers')
        vers = 'v. ' + sys.argv[ind+1]
    else:
        vers = 'unknown version'

    target, base = [], []
    if 'py' in sys.argv or 'all' in sys.argv:
        target.append('py')
        base.append('pybase.py')
    if 'cpp' in sys.argv or 'all' in sys.argv:
        target.append('cpp')
        base.append('cppbase.hpp')
    if 'java' in sys.argv or 'all' in sys.argv:
        target.append('java')
        base.append('javabase.java')
    if 'm' in sys.argv or 'all' in sys.argv:
        target.append('m')
        base.append('mbase.m')
    if 'cs' in sys.argv or 'all' in sys.argv:
        target.append('cs')
        base.append('csbase.cs')

    for t, b in zip(target, base):
        omask = MaskParser(mask)
        generator = commongen.Generator.factory(t)
        generator.parse(omask)
        generator.write_to_file(b, odir, libpath, exepath, vers)
