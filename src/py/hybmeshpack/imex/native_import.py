from hybmeshpack.basic.interf import SilentCallbackCancel2
from hybmeshpack.hmcore import hmxml


def read_contour2(doc, node, cb=None):
    " -> c object, object name"
    r = hmxml.import_contour2(doc, node, cb)
    if r is not None:
        return r[0], r[1]


def read_grid2(doc, node, cb=None):
    " -> c object, object name"
    r = hmxml.import_grid2(doc, node, cb)
    if r is not None:
        return r[0], r[1]


def read_surface3(doc, node, cb=None):
    " -> c object, object name"
    r = hmxml.import_surface3(doc, node, cb)
    if r is not None:
        return r[0], r[1]


def read_grid3(doc, node, cb=None):
    " -> c object, object name"
    r = hmxml.import_grid3(doc, node, cb)
    if r is not None:
        return r[0], r[1]


def contour2(doc, node, name=None, cb=None):
    s = "CONTOUR2D"
    if name:
        s = s + "[@name='%s']" % name
    q = hmxml.query(node, s, required='>0')[0]
    return read_contour2(doc, q, cb)


def grid2(doc, node, name=None, cb=None):
    s = "GRID2D"
    if name:
        s = s + "[@name='%s']" % name
    q = hmxml.query(node, s, required='>0')[0]
    return read_grid2(doc, q, cb)


def surface3(doc, node, name=None, cb=None):
    s = "SURFACE3D"
    if name:
        s = s + "[@name='%s']" % name
    q = hmxml.query(node, s, required='>0')[0]
    return read_surface3(doc, q, cb)


def grid3(doc, node, name=None, cb=None):
    s = "GRID3D"
    if name:
        s = s + "[@name='%s']" % name
    q = hmxml.query(node, s, required='>0')[0]
    return read_grid3(doc, q, cb)


def all_contours2(doc, node, cb=None):
    " -> c2 cdata, c2names"
    if cb is None:
        cb = SilentCallbackCancel2()
    c2, c2n = [], []
    qs = hmxml.query(node, "CONTOUR2D")
    for i, q in enumerate(qs):
        cb1 = cb.subcallback(i, len(qs))
        c, n = read_contour2(doc, q, cb1)
        c2.append(c)
        c2n.append(n)
    return c2, c2n


def all_grids2(doc, node, cb=None):
    " -> g2 cdata, g2names"
    if cb is None:
        cb = SilentCallbackCancel2()
    g2, g2n = [], []
    qs = hmxml.query(node, "GRID2D")
    for i, q in enumerate(qs):
        cb1 = cb.subcallback(i, len(qs))
        g, n = read_grid2(doc, q, cb1)
        g2.append(g)
        g2n.append(n)
    return g2, g2n


def all_surfaces3(doc, node, cb=None):
    " -> s3 cdata, s3names"
    if cb is None:
        cb = SilentCallbackCancel2()
    s3, s3n = [], []
    qs = hmxml.query(node, "SURFACE3D")
    for i, q in enumerate(qs):
        cb1 = cb.subcallback(i, len(qs))
        s, n = read_surface3(doc, q, cb1)
        s3.append(s)
        s3n.append(n)
    return s3, s3n


def all_grids3(doc, node, cb=None):
    " -> g3 cdata, g3names"
    if cb is None:
        cb = SilentCallbackCancel2()
    g3, g3n = [], []
    qs = hmxml.query(node, "GRID3D")
    for i, q in enumerate(qs):
        cb1 = cb.subcallback(i, len(qs))
        g, n = read_grid3(doc, q, cb1)
        g3.append(g)
        g3n.append(n)
    return g3, g3n


def all_geom(doc, node, cb=None):
    """ -> [conts2], [grids2], [surfaces3], [grids3],
           [c2names], [g2names], [s3names], [g3names]
    """
    # contours
    c2, c2n = all_contours2(doc, node, cb)
    # surfaces
    s3, s3n = all_surfaces3(doc, node, cb)
    # grids2
    g2, g2n = all_grids2(doc, node, cb)
    # grids3
    g3, g3n = all_grids3(doc, node, cb)
    return c2, g2, s3, g3, c2n, g2n, s3n, g3n
