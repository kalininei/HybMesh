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