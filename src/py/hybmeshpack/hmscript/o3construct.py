"3D objects construction routines"
from hybmeshpack import com
from hybmeshpack.hmscript import flow, hmscriptfun
import copy
from datachecks import (icheck, UListOr1, Bool, Point2D, ZType, Grid2D, UInt,
                        Float, NoneOr, Func, Grid3D, ASurf3D, Or, IncList)


@hmscriptfun
def grid3_bnd_to_surface(gid, separate=False):
    """ Returns surface object built out of grid boundary

    :param gid: 3D grid identifier

    :param bool separate: whether grid surface should be separated
       into singly connected set of surfaces

    :returns: grid identifier if **separate** is False or
       list of grid identifiers otherwise

    """
    icheck(0, Grid3D())
    icheck(1, Bool())

    c = com.surfcom.Grid3BndToSurface({"grid_name": gid,
                                      "separate": separate})
    flow.exec_command(c)
    if separate:
        return c.added_surfaces3()
    else:
        return c.added_surfaces3()[0]


@hmscriptfun
def tetrahedral_fill(domain):
    """ Fills 3D domain with tetrahedral mesh

    :param domain: surface/3d grid identifier (or list of identifiers)

    :returns: 3d grid identifier

    Domain is defined by any number of closed surfaces passed in **domain**
    argument. Internal nesting procedure
    will be executed to built target domain out of given surfaces.
    Boundary surface could possibly contain faces built
    by any number of vertices.
    However if boundary face is not a triangle than a n-side pyramid
    will be built at its site. Hence resulting grid will not be strictly
    tetrahedral.

    .. note::

        By now program uses simplified bounding box based nesting
        algorithm. It could give improper results for complicated
        surface structures. Be sure that passed surface list nesting equals
        nesting of respective surface bounding boxes.
    """
    icheck(0, UListOr1(ASurf3D()))

    if not isinstance(domain, list):
        domain = [domain]

    c = com.grid3dcom.TetrahedralFill({"source": domain})
    flow.exec_command(c)
    return c.added_grids3()[0]


@hmscriptfun
def extrude_grid(obj, zcoords, bottombc=0, topbc=0, sidebc=None):
    """ Creates 3D grid by extrusion of 2D grid along z-axis

    :param obj: 2d grid identifier

    :param list-of-floats zcoords: increasing vector of z values
      which will be used to create 3d points

    :param bottombc:

    :param topbc: values which define boundary features of
      3d grid at ``z=min(zcoords)`` and ``z=max(zcoords)``
      surfaces respectively.
      Could be either a single boundary identifier for a whole
      surface or a function: ``(float x, float y, int cell_index)->bindex``
      which takes central cell point x, y
      coordinates and cell index as arguments and returns boundary type
      (see example below).

    :param sidebc: defines boundary features for side surfaces.

      * If None than boundary types will be taken from corresponding
        edges of 2D grid
      * If single boundary identifier then whole side surface will
        have same boundary type

    :returns: 3D grid identifier

    Use :func:`partition_segment` to define non-equidistant
    **zcoords** with any desired refinement.

    Example:

      .. literalinclude:: ../../testing/py/fromdoc/ex_extrude.py
          :start-after: vvvvvvvvvvvvvvvvvvvvvvvv
          :end-before: ^^^^^^^^^^^^^^^^^^^^^^^^
    """
    icheck(0, Grid2D())
    icheck(1, IncList(Float()))
    icheck(2, Or(ZType(), Func(nargs=3)))
    icheck(3, Or(ZType(), Func(nargs=3)))
    icheck(4, NoneOr(ZType()))

    # calculate boundary types
    if not isinstance(bottombc, int) or not isinstance(topbc, int):
        grid = flow.receiver.get_grid2(obj)
        cc_pnt = grid.raw_data('centers')
    if isinstance(bottombc, int):
        bbot = [bottombc]
    elif callable(bottombc):
        it = iter(cc_pnt)
        bbot = [bottombc(p[0], p[1], i) for i, p in enumerate(zip(it, it))]
    if isinstance(topbc, int):
        btop = [topbc]
    elif callable(topbc):
        it = iter(cc_pnt)
        btop = [topbc(p[0], p[1], i) for i, p in enumerate(zip(it, it))]
    if sidebc is None:
        bside = None
    elif isinstance(sidebc, int):
        bside = sidebc

    c = com.grid3dcom.ExtrudeZ({"base": obj,
                                "zvals": copy.deepcopy(zcoords),
                                "bside": bside,
                                "btop": btop,
                                "bbot": bbot})
    flow.exec_command(c)
    return c.added_grids3()[0]


@hmscriptfun
def revolve_grid(obj, p1, p2, n_phi=None,
                 phi=None, btype1=0, btype2=0, merge_central=False):
    """ Creates 3D grid by revolution of 2D grid around a vector

    :param obj: 2d grid identifier

    :param p1:

    :param p2: points in [x, y] format which define vector of rotation

    :param int n_phi: partition along circular coordinate.
       If this parameter is defined then [0, 360] range will be divided
       into equal parts and full revolution solid will be build.

    :param list-of-floats phi: increasing vector defining
       custom partition of angular range.
       This parameter will be processed if **n_phi** is None.
       If the last value of **phi** is not equal to first one
       plus 360 degree than
       incomplete revolution solid will be built.

    :param btype1:

    :param btype2: boundary identifiers for surfaces which will be build
       as a result of incomplete rotation at end values of **phi** vector.

    :param bool merge_central: if rotation vector coincides
       with boundary edges of input grid then this parameter
       defines whether central cells derived from the revolution
       of respective boundary cells should be merged into one
       complex finite volume (True) or left as they are (False).

    :returns: 3D grid identifier

    All points of input grid should lie to the one side of rotation
    vector.

    Use :func:`partition_segment` to define non-equidistant
    **phi** with any desired refinement if needed.

    """
    icheck(0, Grid2D())
    icheck(1, Point2D())
    icheck(2, Point2D(noteq=p1))
    icheck(3, NoneOr(UInt(minv=3)))
    icheck(4, NoneOr(IncList(Float())))
    icheck(5, ZType())
    icheck(6, ZType())
    icheck(7, Bool())

    # calculate phi's
    if n_phi is None:
        inp_phi = copy.deepcopy(phi)
    else:
        inp_phi = []
        for i in range(n_phi + 1):
            inp_phi.append(360. * i / n_phi)

    c = com.grid3dcom.Revolve({"base": obj,
                               "p1": p1,
                               "p2": p2,
                               "phi": inp_phi,
                               "bt1": btype1,
                               "bt2": btype2,
                               "center_tri": not merge_central})
    flow.exec_command(c)
    return c.added_grids3()[0]
