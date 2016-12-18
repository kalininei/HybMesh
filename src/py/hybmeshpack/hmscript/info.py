from hybmeshpack.hmscript import flow
from hybmeshpack import hmcore as hmcore
from hybmeshpack.hmcore import libhmcport
from hybmeshpack.hmcore import c2 as c2core
from hybmeshpack.hmcore import g2 as g2core
from hybmeshpack.gdata import grid2
from hybmeshpack import progdata
from hybmeshpack.basic import geom as bg


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
    grid = flow.get_receiver().get_grid(name=g1)[2]
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
       dictionary representing contour length, total number of nodes, edges,
       number of edges in each subcontour and number of edges
       of each boundary type::

         {'Nnodes': int,
          'Nedges': int,
          'subcont': [list-of-int],    # edges number in each subcontour
          'btypes': {btype(int): int}  # boundary type: number of edges
          'length': float
         }

    """
    cont = flow.get_receiver().get_any_contour(c1)
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
    ret['length'] = cont.length()
    return ret


def info_grid3d(g1):
    """ TODO
    """
    g = flow.get_receiver().get_grid3(name=g1)[2]
    ret = {}
    ret['Nnodes'] = g.n_points()
    ret['Nedges'] = g.n_edges()
    ret['Nfaces'] = g.n_faces()
    ret['Ncells'] = g.n_cells()
    return ret


def info_surface(s1):
    """ TODO
    """
    s = flow.get_receiver().get_any_surface(s1)
    ret = {}
    ret['Nnodes'] = s.n_points()
    ret['Nedges'] = s.n_edges()
    ret['Nfaces'] = s.n_faces()
    ret['btypes'] = {}
    for b in s.btypes():
        if b not in ret['btypes']:
            ret['btypes'][b] = 0
        ret['btypes'][b] += 1
    return ret


def get_point(obj, ind=None, vclosest=None, eclosest=None, cclosest=None):
    """TODO
    """
    g = flow.get_receiver().get_any(obj)
    ret = None
    if ind is not None:
        ret = g.points[ind]
    if vclosest is not None:
        p = bg.Point2(vclosest[0], vclosest[1])
        ret = g.points[g.closest_point_index(p)]
    if eclosest is not None:
        p = bg.Point2(eclosest[0], eclosest[1])
        eds = g.edges if isinstance(g, grid2.Grid2) else g.edges_points()
        ret = g.closest_edge_point(p, eds)
    if cclosest is not None:
        p = bg.Point2(cclosest[0], cclosest[1])
        if isinstance(g, grid2.Grid2):
            g = g.cont
        g = g.simplify(0)
        ret = g.closest_edge_point(p, g.edges_points())
    if ret is not None:
        return [ret.x, ret.y]
    else:
        raise ValueError


def registered_contours():
    """Returns list of all contours identifiers
    """
    return flow.get_receiver().get_ucontour_names()


def registered_grids():
    """ Returns list of all grids identifiers
    """
    return flow.get_receiver().get_grid_names()


def registered_grids3d():
    """ Returns list of all 3d grids identifiers
    """
    return flow.get_receiver().get_grid3_names()


def registered_surfaces():
    """Returns list of all surfaces identifiers
    """
    return flow.get_receiver().get_usurface_names()


def registered_btypes():
    """ Returns list of all boundary types as (index, name) tuples
    """
    d = flow.get_receiver().get_bnd_types()._data
    return [(x.index, x.name) for x in d]


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
        cont = flow.get_receiver().get_ucontour(name=c)[2]
    except KeyError:
        cont = flow.get_receiver().get_grid(name=c)[2].cont

    # get c object
    ccont = c2core.cont2_to_c(cont)
    # calculate area
    hmcore.libhmcport.ecollection_area.restype = ct.c_double
    ret = float(libhmcport.ecollection_area(ccont))
    # free c object
    libhmcport.free_ecollection_container(ccont)

    return ret


def pick_contour(pnt, contlist=[]):
    import copy
    """not documented
    """
    if len(contlist) > 0:
        lst = copy.deepcopy(contlist)
    else:
        lst = flow.get_receiver().get_all_names2()
    cpts = [get_point(c, eclosest=pnt) for c in lst]
    cpts = [[x[0] - pnt[0], x[1] - pnt[1]] for x in cpts]
    meas = [x[0] * x[0] + x[1] * x[1] for x in cpts]
    minindex, _ = min(enumerate(meas), key=lambda x: x[1])
    return lst[minindex]


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
    grid = flow.get_receiver().get_grid(name=g)[2]
    skew = g2core.get_skewness(grid, threshold)
    if (len(skew) == 0):
        raise Exception("Can not calculate skewness")
    skew['ok'] = (len(skew['bad_cells']) == 0)
    return skew
