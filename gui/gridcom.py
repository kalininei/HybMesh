import bgeom
import grid2
import command
from unite_grids import unite_grids


class AddUnfRectGrid(command.Command):
    " Add uniform rectangular grid "
    def __init__(self, p0, p1, nx, ny, name):
        if name == "":
            name = "Grid1"
        super(AddUnfRectGrid, self).__init__(p0, p1, nx, ny, name)
        self.p0, self.p1 = p0, p1
        self.Nx, self.Ny = nx, ny
        self.name = name
        self.Grid = None

    @classmethod
    def from_strings(cls, slist):
        p0 = bgeom.Point2(float(slist[0]), float(slist[1]))
        p1 = bgeom.Point2(float(slist[2]), float(slist[3]))
        Nx, Ny = int(slist[4]), int(slist[5])
        name = slist[6]
        return cls(p0, p1, Nx, Ny, name)

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

    def __init__(self, p0, r, na, nr, coef, is_trian, name):
        if name == "":
            name = "CircularGrid1"
        super(AddUnfCircGrid, self).__init__(p0, r, na, nr,
                coef, is_trian, name)
        self.p0 = p0
        self.rad = r
        self.Na, self.Nr = na, nr
        self.coef = coef
        self.is_trian = is_trian
        self.name = name
        self.Grid = None

    @classmethod
    def from_strings(cls, slist):
        p0 = bgeom.Point2(float(slist[0]), float(slist[1]))
        r = float(slist[2])
        Na, Nr = int(slist[3]), int(slist[4])
        coef = float(slist[5])
        is_trian = bool(slist[6])
        name = slist[7]
        return cls(p0, r, Na, Nr, coef, is_trian, name)

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
    def __init__(self, old_name, new_name):
        super(RenameGrid2, self).__init__(old_name, new_name)
        self.oldName, self.newName = old_name, new_name

    #overriden from Command
    @classmethod
    def _method_code(cls):
        return "RenameGrid2"

    @classmethod
    def from_strings(cls, slist):
        return cls(slist[0], slist[1])

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
    def __init__(self, rem_grid):
        ' rem_grid - string name of the removing grid '
        super(RemoveGrid2, self).__init__(rem_grid)
        self.remGrid = rem_grid

    #overriden from Command
    @classmethod
    def _method_code(cls):
        return "RemoveGrid2"

    @classmethod
    def from_strings(cls, slist):
        return cls(slist[0])

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

    def __init__(self, name, buf, den):
        self.name = name
        self.buf = buf
        self.den = den

    def __str__(self):
        return " ".join(map(str, [self.name, self.buf, self.den]))


class UniteGrids(command.Command):
    def __init__(self, grid_name, source_grids):
        ' grid_name - name of the new grid, source_grids - [UniteOpts],  '
        super(UniteGrids, self).__init__(grid_name, *source_grids)
        self.grid_name = grid_name
        self.source = source_grids

    @classmethod
    def _method_code(cls):
        return "UniteGrids"

    @classmethod
    def from_strings(cls, slist):
        it = iter(slist)
        nm = it.next()
        src = []
        for n, b, d in zip(it, it, it):
            src.append(UniteOpts(n, float(b), int(d)))
        return cls(nm, src)

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
        super(MoveGrids, self).__init__(dx, dy, *names)
        self.dx, self.dy = dx, dy
        self.names = names

    @classmethod
    def from_strings(cls, slist):
        [dx, dy] = slist[0:2]
        return cls(float(dx), float(dy), slist[2:])

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

    def __init__(self, x0, y0, angle, names):
        super(RotateGrids, self).__init__(x0, y0, angle, *names)
        self.x0, self.y0, self.angle = x0, y0, angle
        self.names = names

    @classmethod
    def from_strings(cls, slist):
        [x0, y0, angle] = map(float, slist[0:3])
        return cls(x0, y0, angle, slist[3:])

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
