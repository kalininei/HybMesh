import addremove
from hybmeshpack.hmcore import g3 as g3core
from hybmeshpack.gdata.grid3 import Grid3
import comopt as co


class NewGrid3DCommand(addremove.AbstractAddRemove):
    "Command with a new grid addition"
    def __init__(self, argsdict):
        super(NewGrid3DCommand, self).__init__(argsdict)

    def _addrem_grid3(self):
        g = self._build_grid()
        if g is not None:
            return [(Grid3(g), self.get_option('name'))], []
        else:
            return [], []

    #function for overloading
    def _build_grid(self):
        '-> grid3.Grid3.cdata pointer'
        raise NotImplementedError


class ExtrudeZ(NewGrid3DCommand):
    "extrude xy grid in z-direction"
    def __init__(self, argsdict):
        super(ExtrudeZ, self).__init__(argsdict)

    @classmethod
    def _arguments_types(cls):
        """
        btop, bbot are either single value array or
                   int array value for each 2d cell;
        bside is * None if boundary should be taken from grid2;
                 * [single_value] if one boundary type for everything
                 * [ed0_0, ed0_1, ... ed0_n, ed1_0, ed1_1, ..., edm_n]:
                   values for each grid2 edge for each z section
                   starting from bottom to top
        zvals is strictly increasing
        """
        return {'name': co.BasicOption(str),
                'base': co.BasicOption(str),
                'zvals': co.ListOfOptions(co.BasicOption(float)),
                'bbot': co.ListCompressedOption([0]),
                'btop': co.ListCompressedOption([0]),
                'bside': co.NoneOr(co.BasicOption(int), None)
                }

    def _build_grid(self):
        g2 = self.grid_by_name(self.get_option('base'))
        return g3core.extrude(g2.cdata,
                              self.get_option('zvals'),
                              self.get_option('bbot'),
                              self.get_option('btop'),
                              self.get_option('bside'))


class Revolve(NewGrid3DCommand):
    "Revolution of 2d grid"

    def __init__(self, argsdict):
        super(Revolve, self).__init__(argsdict)

    @classmethod
    def _arguments_types(cls):
        return {'name': co.BasicOption(str, None),
                'base': co.BasicOption(str),
                'p1': co.Point2Option(),
                'p2': co.Point2Option(),
                'phi': co.ListOfOptions(co.BasicOption(float)),
                'bt1': co.BasicOption(int, 0),
                'bt2': co.BasicOption(int, 0),
                'center_tri': co.BoolOption(True)
                }

    def _build_grid(self):
        g2 = self.grid_by_name(self.get_option('base'))
        vec = [self.get_option('p1'), self.get_option('p2')]
        return g3core.revolve(g2.cdata, vec, self.get_option('phi'),
                              self.get_option('center_tri'),
                              self.get_option('bt1'),
                              self.get_option('bt2'))


class TetrahedralFill(NewGrid3DCommand):
    def __init__(self, argsdict):
        super(TetrahedralFill, self).__init__(argsdict)

    @classmethod
    def _arguments_types(cls):
        return {'name': co.BasicOption(str, None),
                'source': co.ListOfOptions(co.BasicOption(str)),
                'constr': co.ListOfOptions(co.BasicOption(str), []),
                'pts': co.ListOfOptions(co.Point3Option(), []),
                'pts_size': co.ListOfOptions(co.BasicOption(float), []),
                }

    def _build_grid(self):
        src = map(self.any_surface_by_name, self.get_option('source'))
        constr = map(self.any_surface_by_name, self.get_option('constr'))
        cb = self.ask_for_callback()
        return g3core.tetrahedral_fill(
            [x.cdata for x in src], [x.cdata for x in constr],
            self.get_option('pts'), self.get_option('pts_size'), cb)


class Merge(NewGrid3DCommand):
    def __init__(self, argsdict):
        super(Merge, self).__init__(argsdict)

    @classmethod
    def _arguments_types(cls):
        return {'name': co.BasicOption(str, None),
                'src1': co.BasicOption(str),
                'src2': co.BasicOption(str),
                }

    def _build_grid(self):
        g1 = self.grid3_by_name(self.get_option('src1'))
        g2 = self.grid3_by_name(self.get_option('src2'))
        return g3core.merge(g1.cdata, g2.cdata)
