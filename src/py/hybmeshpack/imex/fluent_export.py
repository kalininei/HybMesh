from hybmeshpack.hmcore import g2 as g2core
from hybmeshpack.hmcore import g3 as g3core


def grid2(fname, grid, btypes=None, per_data=None, cb=None):
    """ btypes: {bindex: bname}
        per_data: [btype_per(int), btype_shadow(int), is_reversed(bool),
                   .....]
    """
    g2core.to_msh(grid.cdata, fname, btypes, per_data, cb)


def grid3(fname, grid, btypes=None, per_data=None, cb=None):
    """ btypes: {bindex: bname}
        per_data: [periodic-0, shadow-0, periodic-point-0, shadow-point-0,
                       periodic-1, shadow-1, periodic-point-1, ...]
    """
    g3core.to_msh(grid.cdata, fname, btypes, per_data, cb)
