import ctypes as ct
from . import libhmcport
from hybmeshpack.basic.geom import Point2
from hybmeshpack.gdata.contour2 import Contour2


def cont_to_c(c):
    """ return c pointer to a contours2.Contour2 object
        works for contours of crossgrid library.
        Will be removed
    """
    #fill points
    c_pnt = (ct.c_double * (2 * c.n_points()))()
    for i, p in enumerate(c.points):
        c_pnt[2 * i] = p.x
        c_pnt[2 * i + 1] = p.y
    #fill edges
    edg = c.edges_points()
    c_edg = (ct.c_int * (2 * c.n_edges()))()
    for i, ce in enumerate(edg):
        c_edg[2 * i] = ce[0]
        c_edg[2 * i + 1] = ce[1]

    return libhmcport.contour_construct(
        ct.c_int(c.n_points()),
        ct.c_int(c.n_edges()), c_pnt, c_edg)


def cont2_to_c(cont):
    """ return c pointer to a contours2.AbstractContour2 objects
    """
    # arguments
    c_npnt = ct.c_int(cont.n_points())
    c_pnts = (ct.c_double * (2 * c_npnt.value))()
    c_nedges = ct.c_int(cont.n_edges())
    for i, p in enumerate(cont.points):
        c_pnts[2 * i] = p.x
        c_pnts[2 * i + 1] = p.y

    c_edges = (ct.c_int * (2 * c_nedges.value))()
    for i, [i0, i1, b] in enumerate(cont.edges_points()):
        c_edges[2 * i] = i0
        c_edges[2 * i + 1] = i1

    return libhmcport.create_ecollection_container(c_npnt, c_pnts,
                                                   c_nedges, c_edges)


def cont2_from_c(c_cont):
    """ returns contours2.Contour2 from c pointer to ecollection """
    if c_cont == 0:
        return None
    c_eds, c_pts = ct.POINTER(ct.c_int)(), ct.POINTER(ct.c_double)()
    npts, neds = ct.c_int(), ct.c_int()
    libhmcport.ecollection_edges_info(c_cont,
                                      ct.byref(npts), ct.byref(neds),
                                      ct.byref(c_pts), ct.byref(c_eds))
    # from c arrays to python lists
    pts = []
    edcon = []
    for i in range(0, 2 * npts.value, 2):
        pts.append(Point2(c_pts[i], c_pts[i + 1]))
    for i in range(0, 2 * neds.value, 2):
        edcon.append([c_eds[i], c_eds[i + 1]])
    libhmcport.free_ecollection_edges_info(c_pts, c_eds)
    ret = Contour2.create_from_point_set(pts, edcon)
    return ret


def free_cont2(cont):
    """ frees pointer to ecollection container"""
    libhmcport.free_ecollection_container(cont)


def clip_domain(cc1, cc2, op, simplify):
    """ cc1, cc2 - c pointers to ecollection containers
    op is {1: union, 2: difference, 3: intersection, 4: xor}
    simplify (bool)
    Returns pointer to resulting container or None
    """
    c_op = ct.c_int(op)
    c_simp = ct.c_int(1 if simplify else 0)
    res = libhmcport.domain_clip(cc1, cc2, c_op, c_simp)
    if (res == 0):
        return None
    else:
        return res
