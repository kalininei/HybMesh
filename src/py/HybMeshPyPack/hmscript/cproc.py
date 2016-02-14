from HybMeshPyPack import com, basic
import HybMeshPyPack.com.contcom
from HybMeshPyPack.basic.geom import Point2
from HybMeshPyPack.hmscript import flow, data


def AddRectCont(p0, p1, bnd):
    """ add rectangular contour
        p0, p1 - corner points as [x0, y0], [x1, y1]
        bnd - boundary type >=0,
            integer or list of 4 integers for each segment.
        returns contour identifier
    """
    b = [0, 0, 0, 0]
    if isinstance(bnd, int):
        b = [bnd, bnd, bnd, bnd]
    elif isinstance(bnd, list):
        b = bnd[0:4]
    c = com.contcom.AddRectCont({"p0": Point2(*p0),
        "p1": Point2(*p1), "bnds": b})
    flow.exec_command(c)
    return c._get_added_names()[1][0]


def GridBndToCont(g1, simplify=True):
    """ Extracts grid boundary to user contour
        g1 - grid identifier
        simplify - if true deletes all non-significant edges
        returns contour identifier
    """
    c = com.contcom.GridBndToContour({"grid_name": g1,
        "simplify": simplify})
    flow.exec_command(c)
    return c._get_added_names()[1][0]


def SimplifyContour(cont, simplify=False, angle=0, separate=False):
    """ Separate and simplify user contour
        cont - source contours ids
        simplify - do simplification: make all segments non-collinear
                   Edges will not be splitted if they have different
                   boundary types.
        angle - minimum acute angle between segments
                after simplification (deg)
        separate - separate contour to ones with not-connected geometry
        returns new contour ids
    """
    c = com.contcom.SimplifyContours(
            {"names": [cont],
                "simplify": simplify,
                "angle": angle,
                "separate": separate})
    flow.exec_command(c)
    return c._get_added_names()[1][0]


def UniteContours(conts):
    """
        Unite contours to single contour with complicated connectivity
        conts - list of contours ids to unite
        returns new contour identifier
    """
    c = com.contcom.UniteContours({"sources": conts})
    flow.exec_command(c)
    return c._get_added_names()[1][0]


def AddBoundaryType(index, name="boundary1"):
    """ Creates boundary type
         index - integer index of boundary
         name - string name of the boundary
         returns boundary identifier
    """
    c = com.contcom.EditBoundaryType({"index": index, "name": name})
    flow.exec_command(c)
    return index


def SetBTypeToContour(cont, btps=None, bfun=None):
    """ Mark user or grid contour segments with btypes.
        cont - contour identifier
        btps - list of btypes identifiers for each segment
            or single identifier for the whole contour
        bfun - function: (x0, y0, x1, y1, bt)->btype identifyer or None
            which returns btype from segment endpoints and old boundary type
    """
    cdata = data.get_any_contour(cont)
    args = {"name": cont}
    if btps is not None:
        if not isinstance(btps, list):
            args[btps] = range(cdata.n_edges())
        else:
            for i, b in enumerate(btps):
                if b not in args:
                    args[b] = []
                args[b].append(i)
    elif bfun is not None:
        ep = cdata.edges_points()
        for i, iep in enumerate(ep):
            p0 = cdata.points[iep[0]]
            p1 = cdata.points[iep[1]]
            r = bfun(p0.x, p0.y, p1.x, p1.y, cdata.edge_bnd(i))
            if r is None:
                continue
            if r not in args:
                args[r] = []
            args[r].append(i)

    c = com.contcom.SetBTypeToContour({"btypes": [args]})
    flow.exec_command(c)


def CreateContour(pnts, bnds=0):
    """ Create contour from sequence of points
        pnts - points list as [[x0, y0], [x1, y1], ... ]. If last point
            cooordinate equals first one then contour is considered closed.
        bnds - list of btype identifiers for each contour segment or
            single identifier for the whole contour.
    """
    pt = []
    for p in pnts:
        pt.append(Point2(*p))
    c = com.contcom.CreateContour({"points": pt,
        "bnds": bnds})
    flow.exec_command(c)
    return c._get_added_names()[1][0]

