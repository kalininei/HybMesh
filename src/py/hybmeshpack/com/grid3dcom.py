import objcom
import command
from hybmeshpack import hmcore as hmcore
from hybmeshpack.hmcore import g2 as g2core
from hybmeshpack.hmcore import g3 as g3core


class NewGrid3DCommand(objcom.AbstractAddRemove):
    "Command with a new grid addition"
    def __init__(self, argsdict):
        if "name" not in argsdict:
            argsdict["name"] = "Grid3D_1"
        super(NewGrid3DCommand, self).__init__(argsdict)

    def _addrem_objects(self):
        #backup info
        g = self._build_grid()
        if g is not None:
            return [], [], [], [], [(self.options['name'], g)], []
        else:
            return [], [], [], [], [], []

    #function for overloading
    def _build_grid(self):
        '-> grid3.Grid3'
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
        return {'name': command.BasicOption(str),
                'base': command.BasicOption(str),
                'zvals': command.ListOfOptions(command.BasicOption(float)),
                'bbot': command.ListCompressedOption(),
                'btop': command.ListCompressedOption(),
                'bside': command.NoneOr(command.BasicOption(int))
                }

    def _build_grid(self):
        so = self.options
        # 0) get grid
        grid = self.grid_by_name(so['base'])
        # 1) grid to c
        c_grid = g2core.grid_to_c(grid)
        # 2) vectors
        c_zvals = hmcore.list_to_c(so['zvals'], float)
        c_bbot = hmcore.list_to_c(so['bbot'], int)
        c_btop = hmcore.list_to_c(so['btop'], int)
        # 3) side boundary conditions
        c_btypes = g2core.boundary_types_to_c(grid)
        bsides = so['bside']
        # 4) call main function
        try:
            c_return = g3core.extrude_call(c_grid, c_btypes,
                                           c_zvals, c_bbot, c_btop,
                                           bsides)
            if c_return is None:
                raise command.ExecutionError("Extrusion failed", self)
            return g3core.grid3_from_c(c_return)
        finally:
            g2core.free_c_grid(c_grid)
            g2core.free_boundary_types(c_btypes)


class Revolve(NewGrid3DCommand):
    "revolution around a defined axis"

    def __init__(self, argsdict):
        super(Revolve, self).__init__(argsdict)

    @classmethod
    def _arguments_types(cls):
        return {'name': command.BasicOption(str),
                'base': command.BasicOption(str),
                'p1': command.Point2Option(),
                'p2': command.Point2Option(),
                'phi': command.ListOfOptions(command.BasicOption(float)),
                'bt1': command.BasicOption(int),
                'bt2': command.BasicOption(int),
                'center_tri': command.BoolOption()
                }

    def _build_grid(self):
        so = self.options
        c_grid, c_btypes = 0, 0
        try:
            # 0) get grid
            grid = self.grid_by_name(so['base'])
            # 1) grid to c
            c_grid = g2core.grid_to_c(grid)
            # 2) vectors
            c_phivals = hmcore.list_to_c(so['phi'], float)
            v = [so['p1'].x, so['p1'].y, so['p2'].x, so['p2'].y]
            c_vector = hmcore.list_to_c(v, float)
            # 3) side boundary conditions
            c_btypes = g2core.boundary_types_to_c(grid)
            # 4) call main function
            c_return = g3core.revolve_call(
                c_grid, c_vector, c_phivals, c_btypes,
                so['bt1'], so['bt2'], so['center_tri'])
            return g3core.grid3_from_c(c_return)
        except Exception as e:
            raise command.ExecutionError("Revolution failed", self, e)
        finally:
            g2core.free_c_grid(c_grid) if c_grid != 0 else None
            g2core.free_boundary_types(c_btypes) if c_btypes != 0 else None


class TetrahedralFill(NewGrid3DCommand):
    def __init__(self, argsdict):
        super(TetrahedralFill, self).__init__(argsdict)

    @classmethod
    def _arguments_types(cls):
        return {'name': command.BasicOption(str),
                'source': command.ListOfOptions(command.BasicOption(str)),
                'pts': command.ListOfOptions(command.BasicOption(float)),
                'pts_size': command.ListOfOptions(command.BasicOption(float)),
                'constr': command.ListOfOptions(command.BasicOption(str))
                }

    def _build_grid(self):
        so = self.options
        try:
            # 0) get source surfaces
            surf = []
            for s in so['source']:
                surf.append(self.any_surface_by_name(s).surface3())
            # 1) get constraint surfaces
            constr = []
            for s in so['constr']:
                constr.append(self.any_surface_by_name(s).surface3())
            # 2) call main function
            c_return = g3core.tetrahedral_fill(
                [x.cdata for x in surf],
                [x.cdata for x in constr],
                so['pts'], so['pts_size'])
            return g3core.grid3_from_c(c_return)
        except Exception as e:
            raise command.ExecutionError("Tetrahedral meshing failed", self, e)
        finally:
            pass
