import ast
import bgeom
import grid2
import command
from unite_grids import unite_grids


class AddUnfRectGrid(command.Command):
    "Add uniform rectangular grid "
    def __init__(self, argsdict):
        super(AddUnfRectGrid, self).__init__(argsdict)
        self.p0, self.p1 = argsdict['p0'], argsdict['p1']
        self.Nx, self.Ny = argsdict['nx'], argsdict['ny']
        if 'name' not in argsdict or argsdict['name'] == '':
            self.name = "Grid1"
        else:
            self.name = argsdict['name']
        self.Grid = None

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

    def _exec(self):
        #backup info
        self.Grid = grid2.UnfRectGrid(self.p0, self.p1, self.Nx, self.Ny)
        #evalution
        self.receiver.grids2[self.name] = self.Grid
        return True

    def _clear(self):
        del self.Grid

    def _undo(self):
        #backup current name of the grid and its index in fw dictionary
        i, k = self.receiver.grids2.get_by_value(self.Grid)
        self.__gridName = k
        self.__gridPosition = i
        #delete Grid from fw. It still presents in self.Grid attribute
        del self.receiver.grids2[k]

    def _redo(self):
        self.receiver.grids2.insert(self.__gridPosition, self.__gridName,
                                    self.Grid)


class AddUnfCircGrid(command.Command):
    "Add uniform circular grid "

    def __init__(self, argsdict):
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
        if 'name' not in argsdict or argsdict['name'] == '':
            argsdict['name'] = 'CircularGrid1'
        self.name = argsdict['name']

        self.Grid = None

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

    def _exec(self):
        #backup info
        self.Grid = grid2.UnfCircGrid(self.p0, self.rad, self.Na, self.Nr,
                self.coef, self.is_trian)
        #evalution
        self.receiver.grids2[self.name] = self.Grid
        return True

    def _clear(self):
        del self.Grid

    def _undo(self):
        #backup current name of the grid and its index in fw dictionary
        i, k = self.receiver.grids2.get_by_value(self.Grid)
        self.__gridName = k
        self.__gridPosition = i
        #delete Grid from fw. It still presents in self.Grid attribute
        del self.receiver.grids2[k]

    def _redo(self):
        self.receiver.grids2.insert(self.__gridPosition, self.__gridName,
                                    self.Grid)


class AddUnfRingGrid(command.Command):
    "Add uniform circular ring"

    def __init__(self, argsdict):
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
        self.name = argsdict['name']

        self.Grid = None

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

    def _exec(self):
        #backup info
        self.Grid = grid2.UnfRingGrid(self.p0, self.irad, self.orad,
                self.Na, self.Nr, self.coef)
        #evalution
        self.receiver.grids2[self.name] = self.Grid
        return True

    def _clear(self):
        del self.Grid

    def _undo(self):
        #backup current name of the grid and its index in fw dictionary
        i, k = self.receiver.grids2.get_by_value(self.Grid)
        self.__gridName = k
        self.__gridPosition = i
        #delete Grid from fw. It still presents in self.Grid attribute
        del self.receiver.grids2[k]

    def _redo(self):
        self.receiver.grids2.insert(self.__gridPosition, self.__gridName,
                                    self.Grid)


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


class RemoveGrid(command.Command):
    def __init__(self, names):
        'name - string name of the removing grid'
        super(RemoveGrid, self).__init__({'names': ' '.join(names)})
        self.names = names
        self.backup = []

    def doc(self):
        return "Remove grids: %s" % ', '.join(self.names)

    #overriden from Command
    @classmethod
    def _method_code(cls):
        return "RemoveGrid"

    @classmethod
    def fromstring(cls, slist):
        a = ast.literal_eval(slist)
        return cls(a['names'].split())

    def _exec(self):
        self.backup = []
        for n in self.names:
            self.backup.append(self.receiver.get_grid(name=n))
            self.receiver.remove_grid(n)
        return True

    def _clear(self):
        self.backup = []

    def _undo(self):
        for v in self.backup:
            self.receiver.grids2.insert(*v)

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


class UniteGrids(command.Command):
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

    def _exec(self):
        #basic grid
        self.unitedGrid, _, _ = self._get_grid(0)
        #unification
        for i in range(1, len(self.source)):
            g, b, d = self._get_grid(i)
            self.unitedGrid = unite_grids(self.unitedGrid, g, b,
                    d, self.fix_bnd)
            if (self.unitedGrid is None):
                return False
        #write result to receiver
        self.receiver.grids2[self.grid_name] = self.unitedGrid
        return True

    def _clear(self):
        del self.unitedGrid

    def _undo(self):
        #backup current name of the grid and its index in fw dictionary
        i, k = self.receiver.grids2.get_by_value(self.unitedGrid)
        self.__gridName = k
        self.__gridPosition = i
        #delete Grid from fw. It still presents in self.Grid attribute
        del self.receiver.grids2[k]

    def _redo(self):
        self.receiver.grids2.insert(self.__gridPosition,
                self.__gridName, self.unitedGrid)


class MoveGrids(command.Command):
    " Move grids "

    def __init__(self, dx, dy, names):
        a = {'dx': dx, 'dy': dy, 'names': ' '.join(names)}
        super(MoveGrids, self).__init__(a)
        self.dx, self.dy = dx, dy
        self.names = names

    def doc(self):
        return "Move grids: " + ", ".join(self.names)

    @classmethod
    def fromstring(cls, slist):
        a = ast.literal_eval(slist)
        return cls(float(a['dx']), float(a['dy']), a['names'].split())

    @classmethod
    def _method_code(cls):
        return "MoveGrids"

    def _exec(self):
        for g in self.names:
            self.receiver.grids2[g].move(self.dx, self.dy)
        return True

    def _clear(self):
        pass

    def _undo(self):
        for g in self.names:
            self.receiver.grids2[g].move(-self.dx, -self.dy)

    def _redo(self):
        self._exec()


class RotateGrids(command.Command):
    " Rotate grids "

    def __init__(self, p0, angle, names):
        a = {'p0': p0, 'angle': angle, 'names': ' '.join(names)}
        super(RotateGrids, self).__init__(a)
        self.x0, self.y0, self.angle = p0.x, p0.y, angle
        self.names = names

    def doc(self):
        return "Rotate grids: " + ", ".join(self.names)

    @classmethod
    def fromstring(cls, slist):
        a = ast.literal_eval(slist)
        return cls(bgeom.Point2.fromstring(a['p0']),
                float(a['angle']), a['names'].split())

    @classmethod
    def _method_code(cls):
        return "RotateGrids"

    def _exec(self):
        for g in self.names:
            self.receiver.grids2[g].rotate(self.x0, self.y0, self.angle)
        return True

    def _clear(self):
        pass

    def _undo(self):
        for g in self.names:
            self.receiver.grids2[g].rotate(self.x0, self.y0, -self.angle)

    def _redo(self):
        self._exec()


class ScaleGrids(command.Command):
    "Scale grids "

    def __init__(self, p0, xpc, ypc, names):
        a = {'p0': p0, 'xpc': xpc, 'ypc': ypc, 'names': ' '.join(names)}
        super(ScaleGrids, self).__init__(a)
        self.p0 = p0
        self.xpc, self.ypc = xpc, ypc
        self.names = names

    def doc(self):
        return "Scale grids: " + ", ".join(self.names)

    @classmethod
    def fromstring(cls, slist):
        a = ast.literal_eval(slist)
        return cls(bgeom.Point2.fromstring(a['p0']),
                float(a['xpc']), float(a['ypc']), a['names'].split())

    @classmethod
    def _method_code(cls):
        return "ScaleGrids"

    def _exec(self):
        for g in self.names:
            self.receiver.grids2[g].scale(self.p0, self.xpc, self.ypc)
        return True

    def _clear(self):
        pass

    def _undo(self):
        for g in self.names:
            self.receiver.grids2[g].unscale(self.p0, self.xpc, self.ypc)

    def _redo(self):
        self._exec()


class CopyGrid(command.Command):
    "Copy grids"

    def __init__(self, srcnames, copynames):
        "srcnames, copynames are the list of grid names"
        a = {'srcnames': ' '.join(srcnames),
                'copynames': ' '.join(copynames)}
        super(CopyGrid, self).__init__(a)
        self.srcnames = srcnames
        self.copynames = copynames
        #new grids as a tuple (ind, name, grid)
        self.created_grids = []

    def doc(self):
        return "Copy grids: " + ", ".join(self.srcnames)

    @classmethod
    def fromstring(cls, slist):
        a = ast.literal_eval(slist)
        return cls(a['srcnames'].split(), a['copynames'].split())

    @classmethod
    def _method_code(cls):
        return "CopyGrid"

    def _exec(self):
        self.created_grids = []
        for name1, name2 in zip(self.srcnames, self.copynames):
            gold = self.receiver.get_grid(name=name1)[2]
            gnew = gold.deepcopy()
            self.receiver.add_grid(name2, gnew)
            self.created_grids.append(self.receiver.get_grid(grid=gnew))
        return True

    def _clear(self):
        self.created_grids = []

    def _undo(self):
        for _, n, _ in self.created_grids:
            self.receiver.remove_grid(n)

    def _redo(self):
        for v in self.created_grids:
            self.receiver.grids2.insert(*v)

