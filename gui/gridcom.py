import ast
import bgeom
import grid2
import command
from unite_grids import unite_grids


class AddUnfRectGrid(command.Command):
    " Add uniform rectangular grid "
    #def __init__(self, p0, p1, nx, ny, name):
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
    " Add uniform circular grid "

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


class RenameGrid2(command.Command):
    def __init__(self, argsdict):
        super(RenameGrid2, self).__init__(argsdict)
        self.oldName, self.newName = argsdict['old'], argsdict['new']

    #overriden from Command
    @classmethod
    def _method_code(cls):
        return "RenameGrid2"

    @classmethod
    def fromstring(cls, slist):
        a = ast.literal_eval(slist)
        return cls(a)

    def _exec(self):
        self.receiver.grids2.changeKey(self.oldName, self.newName)
        return True

    def _clear(self):
        pass

    def _undo(self):
        self.receiver.grids2.changeKey(self.newName, self.oldName)

    def _redo(self):
        self._exec()


class RemoveGrid2(command.Command):
    def __init__(self, name):
        ' name - string name of the removing grid '
        super(RemoveGrid2, self).__init__({'name': name})
        self.remGrid = name

    #overriden from Command
    @classmethod
    def _method_code(cls):
        return "RemoveGrid2"

    @classmethod
    def fromstring(cls, slist):
        a = ast.literal_eval(slist)
        return cls(a['name'])

    def _exec(self):
        self.backupIndex, self.backupGrid = \
            self.receiver.grids2.get_by_key(self.remGrid)
        del self.receiver.grids2[self.remGrid]
        return True

    def _clear(self):
        del self.backupGrid

    def _undo(self):
        self.receiver.grids2.insert(self.backupIndex,
                                    self.remGrid, self.backupGrid)

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
            args[s0] = UniteOpts,
            args[s1] = UniteOpts,
            ....
        """
        super(UniteGrids, self).__init__(kwargs)
        self.grid_name = kwargs['name']
        self.source = [None for i in range(len(kwargs) - 1)]
        for k, v in kwargs.items():
            if k[0] == 's':
                num = int(k[1:])
                self.source[num] = v

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
            self.unitedGrid = unite_grids(self.unitedGrid, g, b, d)
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
