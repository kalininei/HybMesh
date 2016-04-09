from hybmeshpack.hmscript import data
from hybmeshpack import hmcore as hmcore
from hybmeshpack.hmcore import libhmcport
from hybmeshpack.hmcore import c2 as c2core
from hybmeshpack.hmcore import g2 as g2core
from hybmeshpack import progdata


def info_grid(g1):
    """Get grid structure information

    Args:
       g1: grid identifier

    Returns:
       dictionary which represents
       total number of nodes/cells/edges
       and number of cells of each type::

         {'Nnodes': int,
          'Ncells': int,
          'Nedges': int,
          'cell_types': {int: int}  # cell dimension: number of such cells
          }

    """
    grid = data.get_grid(name=g1)[2]
    ret = {}
    ret['Nnodes'] = grid.n_points()
    ret['Ncells'] = grid.n_cells()
    ret['Nedges'] = grid.n_edges()
    ret['cell_types'] = grid.cell_types_info()
    return ret


def info_contour(c1):
    """Get contour structure information

    Args:
       c1: contour or grid identifier

    Returns:
       dictionary representing total number of nodes, edges,
       number of edges in each subcontour and number of edges
       of each boundary type::

         {'Nnodes': int,
          'Nedges': int,
          'subcont': [list-of-int],    # edges number in each subcontour
          'btypes': {btype(int): int}  # boundary type: number of edges
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
    """Returns list of all contours identifiers
    """
    return data.get_ucontour_names()


def registered_grids():
    """ Returns list of all grids identifiers
    """
    return data.get_grid_names()


def registered_btypes():
    """ Returns list of all boundary types as (index, name) tuples
    """
    return [(x.index, x.name) for x in data.get_bnd_types()._data]


def domain_area(c):
    """Calculates area of the domain bounded by the ``c`` contour

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
    ccont = c2core.cont2_to_c(cont)
    # calculate area
    hmcore.libhmcport.ecollection_area.restype = ct.c_double
    ret = float(libhmcport.ecollection_area(ccont))
    # free c object
    libhmcport.free_ecollection_container(ccont)

    return ret


def check_compatibility(vers, policy=1):
    """Checks version compatibility. Notifies if current version
    of hymbesh is not fully compatible with input version.

    Args:
       vers (str): version in ``"0.1.2"`` format

    Kwargs:
       policy (int): if versions are incompatible then:

       * 0 - do nothing
       * 1 - report warning to cout
       * 2 - raise Exception

    Returns:
        False if versions are incompatible, True otherwise.

    """
    versc = progdata.HybMeshVersion.current()
    versi = progdata.HybMeshVersion(vers)
    # last checked for 0.2.1 version
    ret = True
    if (versi > versc):
        ret = False
    if (versi < progdata.HybMeshVersion("0.2.1")):
        ret = False
    if not ret:
        message = "Current version of hybmesh %s is not " \
            "fully compatible with version %s" % (versc, vers)
        if policy == 1:
            print "WARNING:", message
        elif policy == 2:
            raise Exception(message)
    return ret


def skewness(g, threshold=0.7):
    """Reports equiangular skewness of grid cells

    Args:
       g: grid identifier

    Kwargs:
       threshold (float in [0, 1]): cells with skewness higher than this
       value are considered bad and will be reported.
       Set it to -1 to get skewness for each cell.

    Returns:
       dictionary with keys::

          {'ok': bool,                  # True if no bad_cells were found
           'max_skew': float,           # maximum skew value in grid
           'max_skew_cell': int,        # index of cell with maximum skew
           'bad_cells': [list-of-int],  # list of bad cell indices
           'bad_skew': [list-of-float]  # list of bad cell skew values
          }

       `bad_cells` and `bad_skew` lists entries correspond to same bad cells

    """
    grid = data.get_grid(name=g)[2]
    skew = g2core.get_skewness(grid, threshold)
    if (len(skew) == 0):
        raise Exception("Can not calculate skewness")
    skew['ok'] = (len(skew['bad_cells']) == 0)
    return skew
