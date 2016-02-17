from HybMeshPyPack.hmscript import data
from HybMeshPyPack.com import cobj


def info_grid(g1):
    """Get grid connectivity information

    Args:
       g1: grid identifier

    Returns:
       dictionary which represents
       total number of nodes/cells/edges
       and number of cells of each type:

       {'Nnodes': int,

        'Ncells': int,

        'Nedges': int,

        'cells_types': {cells_dims: n1, ... }}

    """
    grid = data.get_grid(name=g1)[2]
    ret = {}
    ret['Nnodes'] = grid.n_points()
    ret['Ncells'] = grid.n_cells()
    ret['Nedges'] = grid.n_edges()
    ret['cells_types'] = grid.cell_types_info()
    return ret


def info_contour(c1):
    """Get contour connectivity information

    Args:
       c1: contour or grid identifier

    Returns:
       dictionary representing total number of nodes, edges,
       number of edges in each subcontour, number of edges of each boundary
       type

       {'Nnodes' : int,

        'Nedges' : int,

        'subcont': [list-of-int],

        'btypes' : {btype(int): int}
       }

    """
    cont = data.get_any_contour(c1)
    ret = {}
    ret['Nnodes'] = cont.n_points()
    ret['Nedges'] = cont.n_edges()
    se = cont.sorted_edges()
    ret['subcont'] = map(len, se)
    ret['btypes'] = {}
    si = [cont.edge_bnd(i) for i in range(cont.n_edges())]
    for s in si:
        if s not in ret['btypes']:
            ret['btypes'][s] = 1
        else:
            ret['btypes'][s] += 1
    return ret


def registered_contours():
    """Get all contours

    Returns:
       List of all contours identifiers

    """
    return data.get_ucontour_names()


def registered_grids():
    """Get all grids

    Returns:
       List of all grids identifiers

    """
    return data.get_grid_names()


def registered_btypes():
    """Get all boundary types

    Returns:
       List of all boundary types as (index, name) tuples

    """
    return [(x.index, x.name) for x in data.get_bnd_types()._data]


def domain_area(c):
    """Calculates area of the domain bounded by contour

    Args:
       c: grid or contour identifier

    Returns:
       positive float or zero for open contours
    """
    import ctypes as ct
    # get contour by name
    try:
        cont = data.get_ucontour(name=c)[2]
    except KeyError:
        cont = data.get_grid(name=c)[2].cont

    # get c object
    ccont = cobj.cont2_to_c(cont)
    # calculate area
    lib_fa = cobj.cport_lib()
    lib_fa.ecollection_area.restype = ct.c_double
    ret = float(lib_fa.ecollection_area(ccont))
    # free c object
    lib_fa.free_ecollection_container(ccont)

    return ret
