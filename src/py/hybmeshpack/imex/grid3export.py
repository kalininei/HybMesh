import hybmeshpack.basic.interf as interf
import hybmeshpack.hmcore as hmcore
import hybmeshpack.hmcore.g3 as g3core


def vtk(grid, fname, cb=None):
    """ cb -- Callback.CB_CANCEL2 callback object or None
    """
    if cb is None:
        cb = interf.SilentCallbackCancel2()
    g3core.to_vtk(grid.cdata, fname, cb)


def vtk_surface(grid, fname, cb=None):
    if cb is None:
        cb = interf.SilentCallbackCancel2()
    g3core.surface_to_vtk(grid.cdata, fname, cb)


def msh(grid, fname, btypes=None, cb=None, per_data=None):
    # callback
    if cb is None:
        cb = interf.SilentCallbackCancel2()

    # boundary names
    if btypes is not None:
        c_bnames = hmcore.boundary_names_to_c(btypes)
    else:
        c_bnames = None

    # periodic
    if per_data is not None:
        if len(per_data) % 4 != 0:
            raise Exception("Invalid length of periodic data array")
        periodic_list = []
        for i in range(len(per_data) / 4):
            [tp1, tp2, p1, p2] = per_data[4 * i:4 * i + 4]
            if not isinstance(tp1, int) or\
                    not isinstance(tp2, int) or\
                    not isinstance(p1, list) or\
                    not isinstance(p2, list) or\
                    len(p1) != 3 or len(p2) != 3:
                raise Exception("Invalid periodic data")
            periodic_list.append(float(tp1))
            periodic_list.append(float(tp2))
            periodic_list.extend(p1)
            periodic_list.extend(p2)
        c_per = hmcore.list_to_c(periodic_list, float)
    else:
        c_per = None

    try:
        g3core.to_msh(grid.cdata, fname, c_bnames, c_per, cb)
    finally:
        if c_bnames is not None:
            hmcore.free_boundary_names(c_bnames)
