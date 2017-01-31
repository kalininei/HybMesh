" 2D object information"
from hybmeshpack.hmscript import flow
from hybmeshpack.hmcore import c2 as c2core
from hybmeshpack.gdata.cont2 import Contour2, closest_contour


def info_grid(gid):
    """Get grid structure information

    :param gid: grid identifier

    :returns: dictionary which represents
       total number of nodes/cells/edges
       and number of cells of each type::

         {'Nnodes': int,
          'Ncells': int,
          'Nedges': int,
          'cell_types': {int: int}  # cell dimension: number of such cells
          }

    """
    grid = flow.receiver.get_grid2(gid)
    ret = {}
    ret['Nnodes'] = grid.n_points()
    ret['Ncells'] = grid.n_cells()
    ret['Nedges'] = grid.n_edges()
    ret['cell_types'] = grid.cell_types_info()
    return ret


def info_contour(cid):
    """Get contour structure information

    :param cid: contour or grid identifier

    :returns:
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
    cont = flow.receiver.get_any_contour(cid).contour2()
    ret = {}
    ret['Nnodes'] = cont.n_points()
    ret['Nedges'] = cont.n_edges()

    sep = [Contour2(s) for s in c2core.quick_separate(cont.cdata)]
    ret['subcont'] = [s.n_edges() for s in sep]

    bt = cont.raw_data('btypes')
    ret['btypes'] = {}
    for s in bt:
        if s not in ret['btypes']:
            ret['btypes'][s] = 1
        else:
            ret['btypes'][s] += 1
    ret['length'] = cont.length()
    return ret


def get_point(obj, ind=None, vclosest=None, eclosest=None, cclosest=None,
              only_contour=True):
    """ Returns object point

    :param obj: grid or contour identifier

    :param int ind: index of point

    :param vclosest:

    :param eclosest:

    :param cclosest: point as [x, y] list

    :param bool only_contour: If that is true then if **objs** is a grid
       then respective grid contour will be used

    :returns: point as [x, y] list

    Only one of **ind**, **vclosest**, **eclosest**, **cclosest**
    arguments should be defined.

    If **ind** is defined then returns point at given index.

    If **vvlosest** point is defined then returns object vertex closest to
    this point.

    If **eclosest** point is defined then returns point owned by an
    object edge closest to input point.

    If **cclosest** point is defined then returns non straight line
    object contour vertex closest to input point.
    """
    try:
        tar = flow.receiver.get_contour2(obj)
    except KeyError:
        tar = flow.receiver.get_grid2(obj)
        if only_contour or cclosest is not None:
            tar = tar.contour()
    if cclosest is not None:
        tar = tar.deepcopy()
        c2core.simplify(tar.cdata, 0, False)
        vclosest = cclosest

    if vclosest is not None:
        return tar.closest_points([vclosest], "vertex")[0]
    elif eclosest is not None:
        return tar.closest_points([eclosest], "edge")[0]
    elif ind is not None:
        return tar.point_at(ind)


def domain_area(cid):
    """Calculates area of the domain bounded by the contour

    :param cid: grid or contour identifier

    :returns: positive float or zero for open contours
    """
    return flow.receiver.get_any_contour(cid).area()


def pick_contour(pnt, contlist=[]):
    """ Returns contour closest to given point

    :param pnt: point as [x, y] list

    :param contlist: list of contour identifier to choose from

    :returns: closest contour identifier

    If **contlist** is empty then looks over all registered contours.
    This procedure does not take 2d grid contours into account.
    """
    conts = map(flow.receiver.get_any_contour, contlist)
    return closest_contour(conts)


def skewness(gid, threshold=0.7):
    """Reports equiangular skewness coefficient (in [0, 1]) of grid cells

    :param gid: grid identifier

    :param float threshold: cells with skewness higher than this
       value are considered bad and will be reported.
       Set it to -1 to get skewness for each cell.

    :return: dictionary with keys::

          {'ok': bool,                  # True if no bad_cells were found
           'max_skew': float,           # maximum skew value in grid
           'max_skew_cell': int,        # index of cell with maximum skew
           'bad_cells': [list-of-int],  # list of bad cell indicies
           'bad_skew': [list-of-float]  # list of bad cell skew values
          }

    `bad_cells` and `bad_skew` lists entries correspond to same bad cells

    """
    skew = flow.receiver.get_grid2(gid).skewness(threshold)
    if (len(skew) == 0):
        raise Exception("Failed to calculate skewness")
    skew['ok'] = (len(skew['bad_cells']) == 0)
    return skew
