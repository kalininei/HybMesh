"""Script for assembling binding headers.
For installation purposes only"""
import sys
import commongen


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


if __name__ == "__main__":
    """possible arguments:
        [py/cpp/java/m/oct/all]
        [odir /path/to/...]
        [libpath /path/to/lib/]
        [exepath /path/to/exe/]
    """
    try:
        dirindex = sys.argv.index('odir')
        odir = sys.argv[dirindex + 1]
    except ValueError:
        odir = 'out'

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
        generator.write_to_file(b, odir, libpath, exepath)
