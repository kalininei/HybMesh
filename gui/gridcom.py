import ast
import copy
import bgeom
import grid2
import command
from unite_grids import unite_grids


class NewGridCommand(command.Command):
    "Command with a new grid addition"
    def __init__(self, argsdict):
        super(NewGridCommand, self).__init__(argsdict)
        self._backup_grid = None
        try:
            self.name = argsdict["name"]
        except KeyError:
            self.name = ""
        if self.name == "":
            self.name = "NewGrid1"

    def _exec(self):
        #backup info
        g = self._build_grid()
        g.build_contour()
        self._backup_grid = self.receiver.add_grid(self.name, g)
        return True

    def _clear(self):
        del self._backup_grid

    def _undo(self):
        self.receiver.remove_grid(self._backup_grid[1])

    def _redo(self):
        self.receiver.grids2.insert(*self._backup_grid)

    #function for overriding
    def _build_grid(self):
        '-> grid2.Grid2'
        raise NotImplementedError


class AddUnfRectGrid(NewGridCommand):
    "Add uniform rectangular grid "
    def __init__(self, argsdict):
        argsdict = copy.deepcopy(argsdict)
        super(AddUnfRectGrid, self).__init__(argsdict)
        self.p0, self.p1 = argsdict['p0'], argsdict['p1']
        self.Nx, self.Ny = argsdict['nx'], argsdict['ny']

    @classmethod
    def fromstring(cls, slist):
        arg = ast.literal_eval(slist)
        arg['p0'] = bgeom.Point2.fromstring(arg['p0'])
        arg['p1'] = bgeom.Point2.fromstring(arg['p1'])
        arg['nx'] = int(arg['nx'])
        arg['ny'] = int(arg['ny'])
        return cls(arg)

    @classmethod
    def _method_code(cls):
        return "AddUnfRectGrid"

    def _build_grid(self):
        return grid2.UnfRectGrid(self.p0, self.p1, self.Nx, self.Ny)


class AddUnfCircGrid(NewGridCommand):
    "Add uniform circular grid "

    def __init__(self, argsdict):
        argsdict = copy.deepcopy(argsdict)
        super(AddUnfCircGrid, self).__init__(argsdict)
        self.p0 = argsdict['p0']
        self.rad = argsdict['rad']
        self.Na, self.Nr = argsdict['na'], argsdict['nr']
        try:
            self.coef = argsdict['coef']
        except KeyError:
            self.coef = 1.0
        try:
            self.is_trian = argsdict['is_trian']
        except KeyError:
            self.is_trian = True

    @classmethod
    def fromstring(cls, slist):
        a = ast.literal_eval(slist)
        a['p0'] = bgeom.Point2.fromstring(a['p0'])
        a['rad'] = float(a['rad'])
        a['na'], a['nr'] = int(a['na']), int(a['nr'])
        if 'coef' in a:
            a['coef'] = float(a['coef'])
        if 'is_trian' in a:
            a['is_trian'] = bool(a['is_trian'])
        return cls(a)

    @classmethod
    def _method_code(cls):
        return "AddUnfCircGrid"

    def _build_grid(self):
        return grid2.UnfCircGrid(self.p0, self.rad, self.Na, self.Nr,
                self.coef, self.is_trian)


class AddUnfRingGrid(NewGridCommand):
    "Add uniform circular ring"

    def __init__(self, argsdict):
        argsdict = copy.deepcopy(argsdict)
        super(AddUnfRingGrid, self).__init__(argsdict)
        self.p0 = argsdict['p0']
        self.irad = argsdict['radinner']
        self.orad = argsdict['radouter']
        self.Na, self.Nr = argsdict['na'], argsdict['nr']
        try:
            self.coef = argsdict['coef']
        except KeyError:
            self.coef = 1.0
        if 'name' not in argsdict or argsdict['name'] == '':
            argsdict['name'] = 'CircularGrid1'

    @classmethod
    def fromstring(cls, slist):
        a = ast.literal_eval(slist)
        a['p0'] = bgeom.Point2.fromstring(a['p0'])
        a['radinner'] = float(a['radinner'])
        a['radouter'] = float(a['radouter'])
        a['na'], a['nr'] = int(a['na']), int(a['nr'])
        if 'coef' in a:
            a['coef'] = float(a['coef'])
        return cls(a)

    @classmethod
    def _method_code(cls):
        return "AddUnfRingGrid"

    def _build_grid(self):
        return grid2.UnfRingGrid(self.p0, self.irad, self.orad,
                self.Na, self.Nr, self.coef)


class RenameGrid(command.Command):
    def __init__(self, oldname, newname):
        a = {"oldname": oldname, "newname": newname}
        super(RenameGrid, self).__init__(a)
        self.oldname, self.newname = oldname, newname
        #actual new grid name can differ from self.newname
        #due to unique names policy. Actual new name is stored
        #in self.backup_newname
        self.backup_newname = None

    #overriden from Command
    def doc(self):
        return "Rename grid: %s" % self.oldname

    @classmethod
    def _method_code(cls):
        return "RenameGrid"

    @classmethod
    def fromstring(cls, slist):
        a = ast.literal_eval(slist)
        return cls(a['oldname'], a['newname'])

    def _exec(self):
        i, _, _ = self.receiver.get_grid(name=self.oldname)
        self.receiver.grids2.change_key(self.oldname, self.newname)
        #backup new name as it can be different from self.newname
        _, self.backup_newname, _ = self.receiver.get_grid(ind=i)
        return True

    def _clear(self):
        pass

    def _undo(self):
        self.receiver.grids2.change_key(self.backup_newname, self.oldname)

    def _redo(self):
        self._exec()


class UniteOpts(object):
    ' Grids unification option: gridname + buffer size + density'

    def __init__(self, name, buf=0, den=5):
        self.name = name
        self.buf = buf
        self.den = den

    def __str__(self):
        d = {'name': self.name,
                'buf': self.buf, 'den': self.den}
        return str(d)

    @classmethod
    def fromstring(cls, slist):
        a = ast.literal_eval(slist)
        return cls(**a)


class UniteGrids(NewGridCommand):
    def __init__(self, **kwargs):
        """ args[name] - name of the new grid,
            args[fix_bnd] - whether to fix all boundary points (=True)
            args[s0] = UniteOpts,
            args[s1] = UniteOpts,
            ....
        """
        super(UniteGrids, self).__init__(kwargs)
        self.grid_name = kwargs['name']
        self.source = [None for i in range(len(kwargs))]
        maxnum = 0
        for k, v in kwargs.items():
            if k[0] == 's':
                num = int(k[1:])
                if num > maxnum:
                    maxnum = num
                self.source[num] = v
        self.source = self.source[:maxnum + 1]
        try:
            self.fix_bnd = kwargs['fix_bnd']
        except KeyError:
            self.fix_bnd = True

    def doc(self):
        ret = "Unite Grids: "
        for s in self.source:
            ret += s.name + " "
        return ret

    @classmethod
    def _method_code(cls):
        return "UniteGrids"

    @classmethod
    def fromstring(cls, slist):
        a = ast.literal_eval(slist)
        for k, v in a.items():
            if k[0] == 's':
                a[k] = UniteOpts.fromstring(v)
        return cls(**a)

    def _get_grid(self, ind):
        """ returns (Grid, buffer size, density) from index """
        g = self.receiver.grids2[self.source[ind].name]
        b = self.source[ind].buf
        d = self.source[ind].den
        return g, b, d

    def _build_grid(self):
        #basic grid
        ret, _, _ = self._get_grid(0)
        #unification
        for i in range(1, len(self.source)):
            g, b, d = self._get_grid(i)
            ret = unite_grids(ret, g, b, d, self.fix_bnd)
            if (ret is None):
                return None
        return ret

